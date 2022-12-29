//=====================================================================
//
//   Mega CD BackupRAM management tool
//
// Copyright (c) 2022 Ian Karlsson
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//=====================================================================

#include <stdio.h>	//printf, snprintf
#include <stdint.h>	//uint8_t, uint16_t, uint32_t
#include <stdlib.h>	//malloc, calloc
#include <string.h>	//memcpy, memcmp
#include <ctype.h>	//isalnum, tolower, toupper
#include <limits.h>	//PATH_MAX

#ifndef NDEBUG
#define DEBUG(...) {fprintf(stdout,__VA_ARGS__);}
#else
#define DEBUG(...) {}
#endif

#define ERROR(...) {fprintf(stderr,__VA_ARGS__); exit(-1);}

struct buram_rs
{
	uint8_t mm; // Bits per symbol
	uint8_t nn; // Symbols per block
	uint8_t index_of[256]; // Antilog lookup
	uint8_t alpha_to[256]; // Log lookup
	const uint8_t *polygen; // Generator polynomial lookup
};

struct buram
{
	struct buram_rs rs8;
	struct buram_rs rs6;

	uint16_t crc_tab[256];

	uint8_t *data;
	uint32_t size;
	uint32_t top_of_dir;

	uint8_t rs_buffer[36]; // Last decoded block
	uint32_t rs_buffer_addr;
	uint16_t rs_status;
};

struct buram_stat
{
	uint8_t format[12];
	uint8_t volume[12];
	uint16_t num_free_blocks;
	uint16_t num_files;
};

struct buram_file
{
	uint8_t filename[11];
	uint8_t use_ecc;
	uint16_t start_block;
	uint16_t size;
};

//=====================================================================
// Initialize Reed-Solomon context
//---------------------------------------------------------------------
static void buram_rs_init(struct buram_rs *rs, uint8_t poly, uint8_t symsize, const uint8_t* polygen)
{
	rs->polygen = polygen;

	rs->mm = symsize;
	rs->nn = (1 << symsize) - 1;

	rs->index_of[0] = 0;
	rs->alpha_to[0] = 0;

	// Generate the Galois field tables
	for(int i = 1, sr = 1; i <= rs->nn; i++)
	{
		rs->index_of[sr] = i;
		rs->alpha_to[i] = sr;
		sr <<= 1;
		if(sr & (1 << symsize))
			sr ^= poly;
		sr &= rs->nn;
	}

}

//=====================================================================
// Add with modulo
//---------------------------------------------------------------------
static uint8_t buram_rs_add_mod(struct buram_rs *rs, int i, uint16_t d)
{
	d += rs->polygen[i];
	while(d >= rs->nn)
		d -= rs->nn;
	return rs->alpha_to[d + 1];
}

//=====================================================================
// Encode Reed-Solomon block
//---------------------------------------------------------------------
static void buram_rs_encode(struct buram_rs *rs, uint8_t *data)
{
	uint8_t *p1 = data + 6;
	uint8_t *p2 = data + 7;
	*p1 = 0;
	*p2 = 0;

	for(int i = 0; i < 6; i++)
	{
		uint8_t d = data[i];
		if(d)
		{
			d = rs->index_of[d] - 1;
			*p1 ^= buram_rs_add_mod(rs, i, d);
			*p2 ^= buram_rs_add_mod(rs, i+6, d);
		}
	}
}

//=====================================================================
// Decode Reed-Solomon block and perform error correction
//---------------------------------------------------------------------
static uint8_t buram_rs_decode(struct buram_rs *rs, uint8_t *data)
{
	uint8_t error_mask = 0;
	uint8_t error_loc = 0;

	// calculate error location (syndrome)
	for(int i = 0; i < 8; i++)
	{
		uint16_t d = data[i];
		error_mask ^= d;
		if(d)
		{
			d = rs->index_of[d] + 6 - i;
			while(d >= rs->nn)
				d -= rs->nn;
			error_loc ^= rs->alpha_to[d + 1];
		}
	}

	// correct a single error
	if(error_mask)
	{
		uint16_t d = rs->nn + rs->index_of[error_loc] - rs->index_of[error_mask];
		while(d >= rs->nn)
			d -= rs->nn;
		if(d < 8)
		{
			data[7 - d] ^= error_mask;
			return 1<<1; // error found and corrected
		}
		return 1<<2; // error found but not corrected
	}
	return 0; // no error found
}


//=====================================================================
// Generate the CRC table using polynomial for CRC-16 CCITT.
//---------------------------------------------------------------------
static void buram_crc_init(struct buram *brm)
{
	for(int i = 0; i < 256; i++)
	{
		uint16_t d = i << 8;
		for(int j = 0; j < 8; j++)
			d = (d << 1) ^ ((d & 0x8000) ? 0x1021 : 0);
		brm->crc_tab[i] = d;
	}
}

//=====================================================================
// Calculate CCITT checksum over 32 bytes
//---------------------------------------------------------------------
static uint16_t buram_crc(struct buram *brm, uint8_t* msg)
{
	uint16_t d = 0;
	for(int i = 0; i < 32; i++)
		d = (d << 8) ^ brm->crc_tab[msg[i] ^ (d >> 8)];
	return d;
}

//=====================================================================
// Interleave 36 bytes of data into 48 bytes
//---------------------------------------------------------------------
static void buram_interleave_data(uint8_t* in, uint8_t* out)
{
	uint16_t sr;
	for(int i = 0; i < 12; i++)
	{
		uint16_t sr = *in++;		// ---- ---- aaaa aaaa
		*out++ = sr;				//           **** **--
		sr <<= 8;
		sr |= *in++;				// aaaa aaaa bbbb bbbb
		*out++ = sr >> 2;			//        ** **** --
		sr <<= 8;
		sr |= *in++;				// bbbb bbbb aaaa aaaa
		*out++ = sr >> 4;			//      **** **--
		*out++ = sr << 2;			//             ** ****
	}
}

//=====================================================================
// Deinterleave 36 bytes of data from 48 bytes
//---------------------------------------------------------------------
static void buram_deinterleave_data(uint8_t* in, uint8_t* out)
{
	uint16_t sr;
	for(int i = 0; i < 12; i++)
	{
		uint16_t sr = *in++;    	// ---- ---- AAAA AAaa
		sr <<= 6;               	// --AA AAAA aa-- ----
		sr = (sr & 0xff00) | *in++;	// --AA AAAA BBBB BBbb
		*out++ = sr >> 6;			//   ** **** **
		sr <<= 6;					// AABB BBBB bb-- ----
		sr = (sr & 0xff00) | *in++;	// AABB BBBB CCCC CCcc
		*out++ = sr >> 4;			//      **** ****
		sr <<= 6;					// BBCC CCCC
		sr = (sr & 0xff00) | *in++;	// BBCC CCCC DDDD DDdd
		*out++ = sr >> 2;			//        ** **** **
	}
}

//=====================================================================
// Interleave the 6-bit RS code
//---------------------------------------------------------------------
// The 6-bit code consists of the most significant 6 bits of each byte in
// the message.  ECC is done on a per-subblock level.
// The order of the bytes in the complete 64-byte block is as follows.
// A* is the first subblock, B* is the second subblock and so on.
// The number indicates the order of the byte within that block.
//
// This diagram shows the assignments of blocks and bytes in the packet.
//       00 01 02 03 04 05 06 07 - BYTE number
//      ------------------------
// 	00 | A0 B0 C0 D0 E0 F0 G0 H0 --- Data bits
// 	08 | D5 A1 B1 C1 D1 E1 F1 G1
// 	10 | H1 E5 A2 B2 C2 D2 E2 F2
// 	18 | G2 H2 F5 A3 B3 C3 D3 E3
// 	20 | F3 G3 H3 G5 A4 B4 C4 D4
// 	28 | E4 F4 G4 H4 H5 A5 B5 C5
//  ---+------------------------
// 	30 | A6 B6 C6 D6 E6 F6 G6 H6 --- Parity bits
// 	38 | A7 B7 C7 D7 E7 F7 G7 H7
static void buram_interleave_rs6(int block, uint8_t *in, uint8_t *out)
{
	static const uint8_t last_data_byte[8] = {0x2d,0x2e,0x2f,0x08,0x11,0x1a,0x23,0x2c};
	for(int i=0; i<5; i++)
		out[block + i*9] = in[i] << 2;
	out[last_data_byte[block]] = in[5] << 2;
	out[0x30 + block] = in[6] << 2;
	out[0x38 + block] = in[7] << 2;
}

//=====================================================================
// Deinterleave the 6-bit RS code so that parity bits are at the bottom
//---------------------------------------------------------------------
static void buram_deinterleave_rs6(int block, uint8_t *in, uint8_t *out)
{
	static const uint8_t last_data_byte[8] = {0x2d,0x2e,0x2f,0x08,0x11,0x1a,0x23,0x2c};
	for(int i=0; i<5; i++)
		out[i] = in[block + i*9] >> 2;
	out[5] = in[last_data_byte[block]] >> 2;
	out[6] = in[0x30 + block] >> 2;
	out[7] = in[0x38 + block] >> 2;
}

//=====================================================================
// Interleave the 8-bit RS code
//---------------------------------------------------------------------
// ECC is done on a per-subblock level. Each block consists of 8 subblocks.
//
// The first byte of the subblock consists of bit 7 of the corresponding
// bytes in the packet. The second byte consists of bit 6, and so on.
// This is done so that bytes 6 and 7 of the subblock (corresponding to
// bits 6 and 7 in the packet) will contain the parity bits.
//
// This diagram shows the assignment of bits within one packet
//       07 06 05 04 03 02 | 01 00 - BIT number
//      -------------------+------
// 	38 | A0 A1 A2 A3 A4 A5 | A6 A7 --- Parity bits
// 	30 | B0 B1 B2 B3 B4 B5 | B6 B7
// 	28 | C0 C1 C2 C3 C4 C5 | C6 C7
// 	20 | D0 D1 D2 D3 D4 D5 | D6 D7
// 	18 | E0 E1 E2 E3 E4 E5 | E6 E7
// 	10 | F0 F1 F2 F3 F4 F5 | F6 F7
// 	08 | G0 G1 G2 G3 G4 G5 | G6 G7
// 	00 | H0 H1 H2 H3 H4 H5 | H6 H7
//              |
//               - Data bits
static void buram_interleave_rs8(int block, uint8_t *in, uint8_t *out)
{
	for(int i=0; i<8; i++)
	{
		uint8_t tmp = in[i];
		for(int j=0; j<8; j++)
		{
			out[block + j*8] = (out[block + j*8] << 1) | (tmp >> 7);
			tmp <<= 1;
		}
	}
}

//=====================================================================
// Deinterleave the 8-bit RS code so that parity bits are at the bottom
//---------------------------------------------------------------------
static void buram_deinterleave_rs8(int block, uint8_t *in, uint8_t *out)
{
	for(int i=0; i<8; i++)
	{
		uint8_t tmp = in[block + i*8];
		for(int j=0; j<8; j++)
		{
			out[j] = (out[j] << 1) | (tmp >> 7);
			tmp <<= 1;
		}
	}
}

//=====================================================================
// Encode an ecc block
//---------------------------------------------------------------------
static void buram_encode(struct buram *brm, uint8_t *buf, uint8_t *out)
{
	uint8_t block_buf[8];

	// Calculate CRC
	uint16_t crc = buram_crc(brm, buf + 2);
	buf[0] = crc >> 8;
	buf[1] = crc;
	buf[34] = ~buf[0];
	buf[35] = ~buf[1];

	// Interleave data
	buram_interleave_data(buf, out);
	// Encode with the 6-bit RS code
	for(int block = 0; block < 8; block++)
	{
		buram_deinterleave_rs6(block, out, block_buf);
		buram_rs_encode(&brm->rs6, block_buf);
		buram_interleave_rs6(block, block_buf, out);
	}
	// Encode with the 8-bit RS code
	for(int block = 0; block < 8; block++)
	{
		buram_deinterleave_rs8(block, out, block_buf);
		buram_rs_encode(&brm->rs8, block_buf);
		buram_interleave_rs8(block, block_buf, out);
	}
}

//=====================================================================
// Decode an ecc block
//---------------------------------------------------------------------
static int buram_decode(struct buram *brm, uint8_t *in, uint8_t *buf)
{
	uint8_t block_buf[8];
	uint16_t flags = 0;

	// Decode the 8-bit RS code
	for(int block = 0; block < 8; block++)
	{
		buram_deinterleave_rs8(block, in, block_buf);
		flags |= buram_rs_decode(&brm->rs8, block_buf);
		buram_interleave_rs8(block, block_buf, in);
	}
	// Decode the 6-bit RS code
	for(int block = 0; block < 8; block++)
	{
		buram_deinterleave_rs6(block, in, block_buf);
		flags |= buram_rs_decode(&brm->rs6, block_buf);
		buram_interleave_rs6(block, block_buf, in);
	}
	// Deinterleave data
	buram_deinterleave_data(in, buf);
	// Calculate CRC
	uint16_t crc = buram_crc(brm, buf + 2);
	uint16_t check_crc1 = (buf[0] << 8) | buf[1];
	uint16_t check_crc2 = ~((buf[34] << 8) | buf[35]);

	if(crc != check_crc1 && crc != check_crc2)
		flags |= 8;
	return flags;
}

//=====================================================================
// Get data from SRAM
//---------------------------------------------------------------------
static int buram_get_data(struct buram *brm, uint32_t address, uint32_t size, uint8_t *buffer)
{
	if(brm->size && brm->data)
		while(size--)
			*buffer++ = brm->data[(address++) % brm->size];
	return size;
}

//=====================================================================
// Get data from SRAM
//---------------------------------------------------------------------
static int buram_set_data(struct buram *brm, uint32_t address, uint32_t size, uint8_t *buffer)
{
	if(brm->size && brm->data)
		while(size--)
			brm->data[(address++) % brm->size] = *buffer++;
	return size;
}

//=====================================================================
// Decode to the temp buffer
//---------------------------------------------------------------------
static uint8_t* buram_decode_buffer(struct buram *brm, uint32_t address)
{
	uint32_t buffer_address = address & -0x40;
	if(brm->rs_buffer_addr != buffer_address)
	{
		uint8_t sram_buf[64];
		buram_get_data(brm, buffer_address, 64, sram_buf);
		brm->rs_status = buram_decode(brm, sram_buf, brm->rs_buffer);
		brm->rs_buffer_addr = buffer_address;
	}
	return brm->rs_buffer + 2 + ((address ^ buffer_address) >> 1);
}

//=====================================================================
// Encode the temp buffer
//---------------------------------------------------------------------
static void buram_encode_buffer(struct buram *brm, uint32_t address, uint8_t *data, uint32_t size)
{
	uint8_t sram_buf[64];
	uint32_t buffer_address = address & -0x40;
	if(size < 32 && brm->rs_buffer_addr != buffer_address)
	{
		buram_get_data(brm, buffer_address, 64, sram_buf);
		brm->rs_status = buram_decode(brm, sram_buf, brm->rs_buffer);
	}
	brm->rs_buffer_addr = buffer_address;
	uint8_t *write_ptr = brm->rs_buffer + 2 + ((address ^ buffer_address) >> 1);
	while(size--)
		*write_ptr++ = *data++;

	buram_encode(brm, brm->rs_buffer, sram_buf);
	buram_set_data(brm, buffer_address, 64, sram_buf);
}

//=====================================================================
// Find a number that occurs at least 3 times
//---------------------------------------------------------------------
static int16_t buram_read_repeat_code(uint8_t* data, int count)
{
	if(data)
	{
		uint16_t buf[count];
		for(int i = 0; i < count; i++)
		{
			buf[i] = (data[0] << 8) | (data[1]);
			data += 2;
		}
		for(int i = 0; i < (count/2); i++)
		{
			int repeats = 0;
			for(int j=i+1; j<count; j++)
			{
				if(buf[i] == buf[j])
					repeats++;
			}
			if(repeats > (count/2))
				return buf[i];
		}
	}
	return -1;
}

//=====================================================================
// Write repeating code
//---------------------------------------------------------------------
static void buram_write_repeat_code(uint8_t* data, int count, uint16_t val)
{
	while(count--)
	{
		*data++ = val >> 8;
		*data++ = val;
	}
}

//=====================================================================
// Get the number of free blocks
//---------------------------------------------------------------------
static int16_t buram_get_free_blocks(struct buram *brm)
{
	if(brm->size >= 64)
		return buram_read_repeat_code(brm->data + brm->top_of_dir + 0x10, 4);
	else
		return -1;
}

//=====================================================================
// Get the number of files in the directory
//---------------------------------------------------------------------
static int16_t buram_get_num_files(struct buram *brm)
{
	if(brm->size >= 64)
		return buram_read_repeat_code(brm->data + brm->top_of_dir + 0x18, 4);
	else
		return -1;
}

//=====================================================================
// Update the free file and block counters
//---------------------------------------------------------------------
static int16_t buram_update_counters(struct buram *brm, int delete, int16_t block_delta)
{
	int16_t num_files = buram_get_num_files(brm);
	int16_t free_blocks = buram_get_free_blocks(brm);
	if(delete)
	{
		num_files--;
		if(num_files & 1)
			free_blocks++;
		free_blocks += block_delta;
	}
	else
	{
		if(num_files & 1)
			free_blocks--;
		num_files++;
		free_blocks -= block_delta;
	}
	uint8_t *ptr = brm->data + brm->top_of_dir + 0x10;
	buram_write_repeat_code(brm->data + brm->top_of_dir + 0x10, 4, free_blocks);
	buram_write_repeat_code(brm->data + brm->top_of_dir + 0x18, 4, num_files);
}

//=====================================================================
// Read a file descriptor
//---------------------------------------------------------------------
static void buram_read_file_desc(uint8_t *buf, struct buram_file *file)
{
	for(int i = 0; i < 11; i++)
		file->filename[i] = buf[i];
	file->use_ecc = buf[11];
	file->start_block = (buf[12] << 8) | buf[13];
	file->size = (buf[14] << 8) | buf[15];
}

//=====================================================================
// Write a file descriptor
//---------------------------------------------------------------------
static void buram_write_file_desc(struct buram_file *file, uint8_t *buf)
{
	for(int i = 0; i < 11; i++)
		buf[i] = file->filename[i];
	buf[11] = file->use_ecc;
	buf[12] = file->start_block >> 8;
	buf[13] = file->start_block;
	buf[14] = file->size >> 8;
	buf[15] = file->size;
}

//=====================================================================
// Initialize backup RAM context
//---------------------------------------------------------------------
int buram_init(struct buram* brm, uint8_t *data, uint32_t size, uint8_t **format_str)
{
	static const uint8_t polygen_rs8[12] = {87, 166, 113, 75, 198, 25, 167, 114, 76, 199, 26, 1};
	static const uint8_t polygen_rs6[12] = {20, 58, 56, 18, 26, 6, 59, 57, 19, 27, 7, 1};
	buram_rs_init(&brm->rs8, 0x1D, 8, polygen_rs8);
	buram_rs_init(&brm->rs6, 3, 6, polygen_rs6);
	buram_crc_init(brm);

	brm->data = data;
	brm->size = 0;
	if(data)
	{
		// Read format
		brm->size = size;
		brm->top_of_dir = size - 0x40;
		// The original MCD BRAM code also checks the format here. We omit that check...
		if(format_str)
		{
			memcpy(*format_str, brm->data + brm->top_of_dir + 0x20, 11);
			*format_str[11] = 0;
		}
		return 0;
	}
	return 1;
}

//=====================================================================
// Return format, number of free blocks and number of files
//---------------------------------------------------------------------
int buram_stat(struct buram* brm, struct buram_stat* stat)
{
	int return_val = 1;
	int16_t free_blocks = -1;
	int16_t num_files = -1;
	if(brm->size)
	{
		int16_t free_blocks = buram_get_free_blocks(brm);
		int16_t num_files = buram_get_num_files(brm);
		if(free_blocks >= 0 && num_files >= 0)
		{
			// Get the format...
			if(stat)
			{
				memcpy(stat->format, brm->data + brm->top_of_dir + 0x20, 11);
				stat->format[11] = 0;
				memcpy(stat->volume, brm->data + brm->top_of_dir, 11);
				stat->volume[11] = 0;
				stat->num_free_blocks = free_blocks;
				stat->num_files = num_files;
			}
			return_val = 0;
		}
	}
	return return_val;
}

//=====================================================================
// Search for a file and return file information
//---------------------------------------------------------------------
int buram_search(struct buram *brm, struct buram_file *file, char *filename)
{
	int files_left = buram_get_num_files(brm);
	uint32_t pos = brm->top_of_dir - 0x20;
	while(files_left > 0)
	{
		uint8_t *buf = buram_decode_buffer(brm, pos);
		if(!memcmp(filename, buf, 11))
		{
			buram_read_file_desc(buf, file);
			return 0;
		}
		pos -= 0x20;
		files_left--;
	}
	return -1;
}

//=====================================================================
// Delete a file
//---------------------------------------------------------------------
int buram_delete(struct buram *brm, char *filename)
{
	struct buram_file file;
	int files_left = buram_get_num_files(brm);
	uint32_t pos = brm->top_of_dir - 0x20;
	while(files_left > 0)
	{
		uint8_t *buf = buram_decode_buffer(brm, pos);
		if(!memcmp(filename, buf, 11))
		{
			buram_read_file_desc(buf, &file);
			uint16_t shift_size = file.size;
			while(files_left > 1)
			{
				uint8_t write_buf[16];
				buf = buram_decode_buffer(brm, pos - 0x20);
				buram_read_file_desc(buf, &file);
				file.start_block -= shift_size;
				buram_write_file_desc(&file, write_buf);
				buram_encode_buffer(brm, pos, write_buf, 0x10);
				memmove(brm->data + file.start_block * 64,
						brm->data + (file.start_block + shift_size) * 64,
						file.size * 64);
				pos -= 0x20;
				files_left--;
			}
			buram_update_counters(brm, 1, shift_size);
			return 0;
		}
		pos -= 0x20;
		files_left--;
	}
	return -1; // No files deleted.
}

//=====================================================================
// Format
//---------------------------------------------------------------------
int buram_format(struct buram *brm, uint8_t *volume)
{
	// The BIOS always sets this to spaces, I think this might be supposed to be the volume name
	// though.
	static const uint8_t initial_volume[11] = "___________";
	// Perhaps this is supposed to indicate the block size, however the BIOS does not use it.
	static const uint8_t initial_unknown[5] = {0x00,0x00,0x00,0x00,0x40};
	// Format and media type
	static const uint8_t initial_format[32] = "SEGA_CD_ROM\x00\x01\x00\x00\x00RAM_CARTRIDGE___";
	if(brm->data)
	{
		// Last block contains the footer, first block is reserved. The remaining block is reserved
		// for the first file in the directory, hence the # of free blocks is the total minus 3.
		uint16_t free_blocks = (brm->size / 64) - 3;
		memcpy(brm->data + brm->top_of_dir, volume ? volume : initial_volume, 11);
		memcpy(brm->data + brm->top_of_dir + 0x0b, initial_unknown, 5);
		buram_write_repeat_code(brm->data + brm->top_of_dir + 0x10, 4, free_blocks);
		buram_write_repeat_code(brm->data + brm->top_of_dir + 0x18, 4, 0);
		memcpy(brm->data + brm->top_of_dir + 0x20, initial_format, 32);
		return 0;
	}
	return -1;
}

//=====================================================================
// Read a file
//---------------------------------------------------------------------
int buram_read(struct buram *brm, uint8_t *output, char *filename, long pos, long size)
{
	struct buram_file file;
	int status = buram_search(brm, &file, filename);
	if(!status)
	{
		if(size == 0)
			size = file.size;
		if(pos + size > file.size)
			size = file.size - pos;

		pos += file.start_block;
		if((pos + size) * 64 > brm->size)
			return -1;

		if(file.use_ecc)
		{
			while(size--)
			{
				uint8_t *buffer = buram_decode_buffer(brm, pos * 64);
				status |= brm->rs_status & 12;
				memcpy(output, buffer, 32);
				output += 32;
				pos++;
			}
		}
		else if(size)
		{
			memcpy(output, brm->data + pos * 64, size * 64);
		}
	}
	return status;
}

//=====================================================================
// Write a file
//---------------------------------------------------------------------
int buram_write(struct buram *brm, struct buram_file *file, uint8_t *input)
{
	// Doing buram_search before buram_delete seems pointless, but it's done
	// so we can check if we have enough blocks before deleting the existing
	// file.
	int free_blocks = buram_get_free_blocks(brm);

	struct buram_file last_file;
	int result = buram_search(brm, &last_file, file->filename);
	if(!result)
	{
		if(free_blocks + last_file.size >= file->size)
			buram_delete(brm, last_file.filename);
		else
			return -1;
	}
	else
	{
		if(free_blocks < file->size)
			return -1;
	}

	int num_files = buram_get_num_files(brm);
	uint32_t pos = brm->top_of_dir - 0x20 - (num_files * 0x20);
	file->start_block = 1; // we start writing from block 1

	if(num_files)
	{
		uint8_t *buf = buram_decode_buffer(brm, pos + 0x20);
		buram_read_file_desc(buf, &last_file);
		file->start_block = last_file.start_block + last_file.size;
	}

	uint8_t write_buf[16];
	buram_write_file_desc(file, write_buf);
	buram_encode_buffer(brm, pos, write_buf, 0x10);

	if(file->use_ecc)
	{
		for(int i = 0; i < file->size; i++)
			buram_encode_buffer(brm, (file->start_block + i) * 64, input + i * 32, 32);
	}
	else
	{
		memcpy(brm->data + file->start_block * 64, input, file->size * 64);
	}

	buram_update_counters(brm, 0, file->size);
	return 0;
}

//=====================================================================
// List directory. Supports wildcard patterns using '*'.
//---------------------------------------------------------------------
int buram_dir(struct buram *brm, struct buram_file *file, char *pattern, long skip, long size)
{
	int size_left = size;
	int files_left = buram_get_num_files(brm);
	uint32_t pos = brm->top_of_dir - 0x20;
	while(files_left > 0)
	{
		uint8_t *buf = buram_decode_buffer(brm, pos);
		int match = 0;
		while(match < 11)
		{
			if(pattern[match] == '*')
				match = 11;
			else if(pattern[match] == buf[match])
				match++;
			else
				break;
		}
		if(match == 11)
		{
			if(skip)
			{
				skip--;
			}
			else if(size_left)
			{
				for(int i = 0; i < 11; i++)
					file->filename[i] = buf[i];
				file->use_ecc = buf[11];
				file->start_block = (buf[12] << 8) | buf[13];
				file->size = (buf[14] << 8) | buf[15];
				file++;
				size_left--;
			}
		}
		pos -= 0x20;
		files_left--;
	}
	return size - size_left;
}

//=====================================================================
// Convert filename from ASCII to BURAM format
//---------------------------------------------------------------------
void buram_ascii_to_text(char* in, uint8_t* out, int size)
{
	for(int i = 0; i < size; i++)
	{
		char c = *in;
		if(c && c != '.')
			in++;
		if((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9'))
			*out++ = c;
		else if(c >= 'a' && c <= 'z')
			*out++ = c - 'a' + 'A';
		else if(c == '*' || c == '?') // Can use '?' to prevent automatic globbing
			*out++ = '*';
		else
			*out++ = '_';
	}
	*out++ = 0;
}

//=====================================================================
// Convert filename from BURAM format to ASCII
//---------------------------------------------------------------------
void buram_text_to_ascii(char* in, uint8_t* out, int size)
{
	// Convert (actually since it's ASCII nothing needs to be done here apart from sanity checking)
	for(int i = 0; i < size; i++)
	{
		char c = *in;
		if(c)
			in++;
		if((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '*')
			*out++ = c;
		else
			*out++ = '_';
	}
	*out = 0;
	// Remove trailing _
	for(int i = size - 1; i > 0; i--)
	{
		if(*(--out) == '_')
			*out = 0;
		else
			break;
	}
}

//=====================================================================
// Open sram image
//---------------------------------------------------------------------
int open(struct buram *brm, const char *file_name, long *file_size, int *file_init, char **file_data)
{
	int status = -1;
	if(*file_init)
	{
		FILE *file = NULL;
		uint8_t *data;

		if(file_name)
			file = fopen(file_name, "rb");

		if(file)
		{
			fseek(file,0,SEEK_END);
			*file_size = ftell(file);
			rewind(file);
			*file_data = malloc(*file_size);
			if(*file_data)
			{
				fread(*file_data, 1, *file_size, file);
				fclose(file);
				*file_init = buram_init(brm, *file_data, *file_size, NULL);
				DEBUG("buram_init() returned %d - loading from '%s'\n", *file_init, file_name);
				status = 0;
			}
		}
		else
		{
			*file_data = calloc(1, *file_size);
			if(*file_data)
			{
				*file_init = buram_init(brm, *file_data, *file_size, NULL);
				DEBUG("buram_init() returned %d - initialized empty file\n", *file_init);
				status = 1;
			}
		}
	}
	return status;
}

//=====================================================================
// Save sram image
//---------------------------------------------------------------------
int save(struct buram *brm, const char *file_name, long file_size, int file_init, char *file_data)
{
	int status = -1;
	if(!file_init && file_data)
	{
		FILE *fd;
		uint8_t *data;

		DEBUG("saving to '%s'...\n", file_name);
		fd = fopen(file_name, "wb"); // Temp.
		if(fd)
		{
			fwrite(file_data, file_size, 1, fd);
			fclose(fd);
			status = 0;
		}
		else
		{
			fprintf(stderr, "failed to open file '%s' for writing\n", file_name);
		}
	}
	return status;
}

//=====================================================================
// Filesystem statistics
//---------------------------------------------------------------------
int stats(struct buram *brm, int verbose)
{
	struct buram_stat stat;
	int status = buram_stat(brm, &stat);

	DEBUG("buram_stat() returned %d\n", status);
	if(verbose)
	{
		char volume[12] = "";
		char format[12] = "";
		char media_id[17] = "";

		buram_text_to_ascii(stat.format, format, 11);
		buram_text_to_ascii(stat.volume, volume, 11);

		if(brm->data && brm->size > 16)
			buram_text_to_ascii(brm->data + brm->size - 16, media_id, 16);

		printf("volume: %s\n",volume);
		printf("format: %s\n",format);
		printf("media_id: %s\n",media_id);
	}
	printf("\nfiles: %d, free blocks: %d\n", stat.num_files, stat.num_free_blocks);
	return status;
}


//=====================================================================
// Run callback function for files matching a pattern
//---------------------------------------------------------------------
int list(struct buram *brm, char* pattern, int fail_if_empty, int verbose,
		int (*callback)(struct buram *brm, struct buram_file *file, int verbose, void *param), void *param)
{
	static const int buf_size = 4;
	struct buram_file files[buf_size];
	int start_pos = 0;
	int num_files;
	int status = 0;
	do
	{
		num_files = buram_dir(brm, files, pattern, start_pos, buf_size);
		if(fail_if_empty && !start_pos && !num_files)
		{
			fprintf(stderr, "no files found\n");
			status |= -1;
		}
		for(int i = 0; i < num_files; i++)
		{
			status |= callback(brm, &files[i], verbose, param);
		}
		start_pos += buf_size;
	}
	while(num_files == buf_size);
	return status;
}

// Callback function to print the directory listing
int cb_list(struct buram *brm, struct buram_file *file, int verbose, void *param)
{
	char ascii_buf[12];
	buram_text_to_ascii(file->filename, ascii_buf, 11);
	if(verbose)
		printf("%c %5d %5d %s\n", file->use_ecc ? 'P' : ' ',
				file->start_block, file->size, ascii_buf);
	else
		printf("%s\n", ascii_buf);
	return 0;
}

// Callback function to extract files
int cb_extract(struct buram *brm, struct buram_file *file, int verbose, void *param)
{
	int status = -1;
	char ascii_buf[12];
	buram_text_to_ascii(file->filename, ascii_buf, 11);
	if(verbose)
		printf("%s\n", ascii_buf);

	uint32_t block_size = file->use_ecc ? 32 : 64;
	uint8_t *buffer = malloc(block_size * file->size);
	if(buffer)
	{
		char filename[PATH_MAX];
		snprintf(filename, PATH_MAX, "%s.bin", ascii_buf);
		int read_status = buram_read(brm, buffer, file->filename, 0, 0);
		if(!read_status)
		{
			printf("opening file '%s' for writing ...\n", filename);
			FILE *fd = fopen(filename, "wb");
			if(fd)
			{
				fwrite(buffer, block_size, file->size, fd);
				fclose(fd);
				status = 0;
			}
			else
			{
				fprintf(stderr, "failed to open file '%s' for writing\n", filename);
			}
			free(buffer);
		}
		else
		{
			fprintf(stderr, "buram_read() failed\n");
		}
	}
	else
	{
		fprintf(stderr, "failed to allocate buffer of size %d\n", block_size * file->size);
	}
	return status;
}

// Callback function to delete files
int cb_delete(struct buram *brm, struct buram_file *file, int verbose, void *param)
{
	char ascii_buf[12];
	buram_text_to_ascii(file->filename, ascii_buf, 11);
	if(verbose)
		printf("%s\n", ascii_buf);
	return buram_delete(brm, file->filename);
}

//=====================================================================
// Add file
//---------------------------------------------------------------------
int add(struct buram *brm, char *path, uint8_t *name, int verbose, int protect)
{
	FILE *fd;
	long size;
	uint8_t *data;
	int status = -1;

	fd = fopen(path, "rb");
	if(fd)
	{
		fseek(fd,0,SEEK_END);
		size = ftell(fd);
		rewind(fd);
		data = calloc(1, size + 63);
		if(data)
		{
			int block_size = protect ? 32 : 64;
			struct buram_file file;
			fread(data, 1, size, fd);
			memcpy(&file.filename, name, 11);
			file.use_ecc = protect ? 0xff : 0;
			file.size = ((size + block_size - 1) & -block_size) / block_size;
			file.start_block = 0;
			status = buram_write(brm, &file, data);
			if (status)
			{
				fprintf(stderr, "%s: buram_write() failed. (%d blocks required, %d available)\n", path, file.size, buram_get_free_blocks(brm));
			}

			free(data);
		}
		else
		{
			fprintf(stderr, "failed to allocate buffer\n");
		}
		fclose(fd);
	}
	else
	{
		fprintf(stderr, "failed to open file '%s' for reading\n", path);
	}
	return status;
}

int main(int argc, char **argv)
{
	struct buram brm;

	int status = 0;
	int dirty = 0;

	int statistics = 0;
	int verbose = 0;
	int protect = 0;
	const char *file_name = "test.brm"; // File used for reading and writing (if write_file is null)
	const char *write_file = 0; // File used for writing

	long file_size = 0x2000;
	int file_init = -1;
	char *file_data;

	char name_buf[12];
	char ascii_buf[12];

	int param, operation;

	for(param = 1, operation = 0; param < argc; param++)
	{
		if(argv[param][0] == '-')
		{
			int this_param, sub_param;
			for(this_param = param, sub_param = 1; argv[this_param][sub_param]; sub_param++)
			{
				switch(argv[this_param][sub_param])
				{
					case 'f': // Set filename
						if(++param < argc)
							file_name = argv[param];
						else
							ERROR("missing parameter for '-f'\n");
						break;
					case 'o': // Set output filename
						if(++param < argc)
							write_file = argv[param];
						else
							ERROR("missing parameter for '-o'\n");
						break;
					case 's': // Set size (untested)
						if(++param < argc)
							file_size = strtoul(argv[param], NULL, 0);
						else
							ERROR("missing parameter for '-s'\n");
						break;
					case 'v': // Verbose mode
						verbose = 1;
						break;
					case 'i': // Print statistics
						statistics = 1;
						break;
					case 'p': // Use protection
						protect = 1;
						break;
					case 'n': // Don't use protection
						protect = 0;
						break;
					default:
						operation = tolower(argv[this_param][sub_param]);
						break;
				}
			}
		}
		else
		{
			if(file_init)
			{
				int open_result = open(&brm, operation == 'c' ? NULL : file_name, &file_size, &file_init, &file_data);
				if((operation == 'c' && open_result < 0) || (operation != 'c' && open_result != 0))
					ERROR("failed to open '%s' for reading (operation=%c, result=%d)\n", file_name, operation, open_result);
			}

			char *fn_start = strrchr(argv[param], '/');
			if(!fn_start)
				strrchr(argv[param], '\\');
			fn_start = fn_start ? fn_start + 1 : argv[param];
			buram_ascii_to_text(fn_start, name_buf, 11);
			buram_text_to_ascii(name_buf, ascii_buf, 11);

			operation = toupper(operation);
			switch(operation)
			{
				case 'C':
					DEBUG("formatting SRAM image\n");
					buram_format(&brm, NULL);
					operation = 'A'; //fall through...
				case 'A':
					DEBUG("adding file '%s' as '%s'...\n",argv[param],ascii_buf);
					status |= add(&brm, argv[param], name_buf, verbose, protect);
					dirty = 1;
					break;
				case 'D':
					DEBUG("deleting file '%s'...\n",ascii_buf);
					status |= list(&brm, name_buf, 0, verbose, cb_delete, NULL);
					dirty = 1;
					break;
				case 'T': // Directory listing
					DEBUG("listing '%s'...\n",ascii_buf);
					status |= list(&brm, name_buf, 0, verbose, cb_list, NULL);
					break;
				case 'X':
					DEBUG("extracting '%s' as '%s.bin'...\n",ascii_buf,ascii_buf);
					status |= list(&brm, name_buf, 1, verbose, cb_extract, NULL);
					break;
			}
		}
	}
	switch(operation)
	{
		case 0:
			if(statistics)
				break;
			// Fall through.
		case 'h':
			printf("buram - Mega CD backup RAM management tool\n");
			printf("version 1.0 (C) 2022 Ian Karlsson\n\n");
			printf("operations:\n");
			printf("\t-c: format new BRAM image and add files\n");
			printf("\t-a [files...]: add files\n");
			printf("\t-d [pattern]: delete files\n");
			printf("\t-t <pattern>: list files\n");
			printf("\t-x <pattern>: extract files\n");
			printf("options:\n");
			printf("\t-f [filename]: set filename for reading and writing\n");
			printf("\t-o [filename]: set filename for writing only\n");
			printf("\t-v: verbose output\n");
			printf("\t-i: display volume information\n");
			printf("\t-s [size]: set size (in bytes) when creating new BRAM image\n");
			printf("\t-p: add files with ECC protection\n");
			printf("\t-n: don't add files with protection\n");
			break;
		case 'c': // Create empty SRAM image
			printf("formatting empty SRAM image\n");
			if(open(&brm, file_name, &file_size, &file_init, &file_data) < 0)
				ERROR("failed to open '%s' for reading\n", file_name);
			buram_format(&brm, NULL);
			dirty = 1;
			break;
		case 't': // List all files
			DEBUG("listing all files...\n");
			if(open(&brm, file_name, &file_size, &file_init, &file_data) != 0)
				ERROR("failed to open '%s' for reading\n", file_name);
			buram_ascii_to_text("*", name_buf, 11);
			list(&brm, name_buf, 0, verbose, cb_list, NULL);
			break;
		case 'x': // Extract all files
			DEBUG("extracting all files...\n");
			if(open(&brm, file_name, &file_size, &file_init, &file_data) != 0)
				ERROR("failed to open '%s' for reading\n", file_name);
			buram_ascii_to_text("*", name_buf, 11);
			list(&brm, name_buf, 1, verbose, cb_extract, NULL);
			break;
		case 'd': // Delete
		case 'a': // Add files
			fprintf(stderr, "no files specified\n");
			break;
	}

	if(statistics)
	{
		if(!open(&brm, file_name, &file_size, &file_init, &file_data) < 0)
			ERROR("failed to open '%s' for reading\n", file_name);
		status |= stats(&brm, verbose);
	}

	if(dirty)
	{
		status |= save(&brm, write_file ? write_file : file_name, file_size, file_init, file_data);
	}

	return status;
}
