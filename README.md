# buram
This program allows you to manipulate Mega CD/Sega CD backup RAM images
using a tar-like interface.

## Compiling
`gcc -o buram buram.c`

## Usage
`<>` indicates optional arguments and `[]` indicates mandatory arguments.

`buram [-operations] <-options> <arguments>`

### Operations
- `-c <files...>`: Format new BRAM image and add files
- `-a [files...]`: Add files
- `-d [pattern]`: Delete files
- `-t <pattern>`: List files
- `-x <pattern>`: Extract files

### Options
- `-f [filename]`: Set filename for reading (and writing if `-o` argument not given)
- `-o [filename]`: Set filename for writing only
- `-v`: Verbose output
- `-i`: Display volume information
- `-s [size]`: Set size in bytes when creating a new BRAM image with `-c`
- `-p` Add files with ECC protection
- `-n` Add files without ECC protection (default)

If no filename is given with `-f` or `-o` options, the default is `test.brm`

### Examples
Create a new archive `myarchive.brm` containing a savefile called `SONICCD`:

	$ buram -cf myarchive.brm SONICCD.bin

Check the number of free blocks in `myarchive.brm`:

	$ buram -if myarchive.brm

	files: 1, free blocks: 114

List all files in `myarchive.brm`:

	$ buram -tf myarchive.brm
	      1    11 SONICCD

The last two operations could be combined like this, additionally adding the `-v` flag to see more information:

	$ buram -tvfi myarchive.brm
		  1    11 SONICCD
	volume: _
	format: SEGA_CD_ROM
	media_id: RAM_CARTRIDGE

	files: 1, free blocks: 114

Extract `SONICCD` from `myarchive.brm`:

	$ buram -xf myarchive.brm SONICCD

Delete `SONICCD` from `myarchive.brm`:

	$ buram -df myarchive.brm SONICCD

### License
&copy; 2022 Ian Karlsson

Licensed under the MIT license, see COPYING

