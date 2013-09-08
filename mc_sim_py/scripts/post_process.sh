mv $1*.log $2/
mv *.*.*$1 $2/
mv *.*.po$1 $2/

$HOME/bin/notify.sh $1
echo "JOB $JOB_NAME FINISHED $1 `date`" >> $HOME/job_finished
