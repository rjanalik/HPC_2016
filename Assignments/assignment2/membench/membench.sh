#!/bin/sh
./membench | sed -e '/:/	s//: /g' -e '/  */	s//	/g' | cut -f2,4,6 > $1.xxx && gnuplot $1.gp
