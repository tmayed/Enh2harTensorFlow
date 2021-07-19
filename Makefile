.PHONY: default
.PHONY: clean

F2PY = f2py3
COMP = --fcompiler=gfortran
FLAGS = --f90flags='-ffree-line-length-1024'
EXT_SUFFIX_LOCAL = .cpython-38-x86_64-linux-gnu.so
NAME = enh_solvers

default: enh

enh: ${NAME}.f90
	${F2PY} ${COMP} ${FLAGS} -m ${NAME} -c $<
	mv ${NAME}$(EXT_SUFFIX_LOCAL) ${NAME}.so

clean :
	rm ${NAME}.so
