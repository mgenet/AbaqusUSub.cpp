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




































