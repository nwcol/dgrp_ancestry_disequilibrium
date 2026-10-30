
log = log/$(Cluster).log
error = log/$(Cluster).$(Process).err
output = log/$(Cluster).$(Process).out
requirements = (Target.HasCHTCStaging == true)
container_image = file:///staging/n/nwcollier/h2py-container.sif

request_cpus   = 1
request_memory = 2GB
request_disk   = 4GB

executable = ex_fet_pvals.sh
arguments = -i1 $(data1) -i2 $(data2) -b $(binsout) -s $(sitesout) --threshold 1e-5

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = \
    compute_block_fet_pvals.py, \
    data/$(data1), \
    data/$(data2)

queue data1,data2,binsout,sitesout from variables.txt
