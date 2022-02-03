import subprocess

prog_name = './lu_decomposition'

for k in range(1, 10):
    for j in range(1, 5):
        for i in range(2, 23, 4):
            command = 'mpiexec -n ' + str(j) + ' ' + prog_name + ' ' + str(i * 100)
            subprocess.call(command, shell=True)
