We will use benchexec to run our experiments: https://github.com/sosy-lab/benchexec.

```bash
# Install benchexec using pip or clone from the repo
pip install benchexec

# Revert to cgroups v1 in case you have ubuntu 2204. If you have lower versions then you can ignore the following 3 commands
echo 'GRUB_CMDLINE_LINUX_DEFAULT="${GRUB_CMDLINE_LINUX_DEFAULT} systemd.unified_cgroup_hierarchy=0"' | sudo tee /etc/default/grub.d/cgroupsv1-for-benchexec.cfg
sudo update-grub

sudo reboot

# Testing that things works
# Note
# 1. replace benchexec with <BENCHEXEC_DIR>/bin/benchexec in case you have cloned it
# 2. Keep a note of capital N. This is the number of concurrent processes used. Small is for name.
PYTHONPATH=scripts/tool-info-modules benchexec -N 4 --read-only-dir / test-data/benchmark-definitions/storm.xml --tool-dir storm/build

PYTHONPATH=scripts/tool-info-modules benchexec -N 4 --read-only-dir / test-data/benchmark-definitions/modest.xml --tool-dir Modest/
```

## Four important things
1. Go through the xml files in the directory `test-data/benchmark-definitions`. These are the input files to the tool benchexec. One needs to specify the tasks to run. Also, one can set the number of required CPUs and memory restriction in this file.
2. Go through the task files in `test-data/tasks`. We will need a script to create these for our need.
3. Go through the results created in the `results` directory. This directory contains the logs as well as the xmls that contains CPU time, wall time, memory measurements, as well as the values we extract from the logs.
4. Go through the tool info modules in the directory `scripts/tool-info-modules`. These are used to create command line as well as parse the log produced by the tool.