import os
os.system("cd lib/ && rm *.so")
os.system("cd lib/ && rm *.o")
compArg_C = ["-L/home/james/gsl/lib", "-I/home/james/gsl/include", "-lgsl", "-lgslcblas", "-lm"]
cArg = ""
if len(compArg_C) > 0 :
    cArg = reduce(lambda v1, v2 : v1+" "+v2, compArg_C)
# ==================================================== #
os.system("cd lib/ && gcc -std=gnu99 -c -fPIC rcdcop.c "+cArg)
os.system("cd lib/ && gcc -std=gnu99 rcdcop.o -shared -o rcdcop.so "+cArg)

