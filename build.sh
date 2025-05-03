#!/bin/bash

gcc -Wall -Wextra \
	-ggdb		  \
	./*.c -I./ 	  \
	-lX11		  \
	-lm			  \
	-o simple_renderer
