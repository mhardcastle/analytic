cd /beegfs/car/mjh/inference-z1.25
export PYTHONPATH=/home/mjh/git/analytic:$PYTHONPATH
python /home/mjh/git/analytic/inference/run_distribution_qsub.py $1 >>run-${1}-out.txt 2>&1
