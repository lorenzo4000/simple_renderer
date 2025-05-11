#!/bin/bash

CFLAGS="-Wall -Wextra -ggdb"

gcc $CFLAGS		  \
	./*.c -I./ 	  \
	-lX11		  \
	-lm			  \
	-o simple_renderer
