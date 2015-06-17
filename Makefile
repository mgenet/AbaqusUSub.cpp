########################################################################
###                                                                  ###
### Created by Martin Genet, 2008-2015                               ###
###                                                                  ###
### Laboratoire de Mécanique et de Technologie (LMT), Cachan, France ###
### Lawrence Berkeley National Laboratory, California, USA           ###
### University of California at San Francisco, USA                   ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

space_ndim = 3

lmt_dir = ../LMT

########################################################################

do: generate
do: compile

########################################################################

generate:
	python usub.cpp.py

compile:
	gcc -c \
            -O3 \
            -std=c++11 \
            -fPIC \
            -DSPACE_NDIM=$(space_ndim) \
            -I. \
            -I$(lmt_dir)/include \
            -o usub.o \
            usub.cpp




































