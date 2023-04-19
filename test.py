#! /usr/bin/env python
import os, subprocess
import shlex
import pprint
# # cmd = shlex.split("env -i bash -c 'source /home/artem/cp2k/tools/toolchain/install/setup'")
# cmd = shlex.split("bash -c 'source /home/artem/cp2k/tools/toolchain/install/setup'")
# proc = subprocess.Popen(cmd, stdout = subprocess.PIPE)
# for line in proc.stdout:
#     # print(line)
#     (key, _, value) = line.partition("=")
#     os.environ[key] = value
# #   print(f'{key} ::: {_} ::: {value}')
# proc.communicate()


pprint.pprint(dict(os.environ))