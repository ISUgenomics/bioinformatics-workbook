SLURM scheduler uses `sbatch` command to submit the jobs. You can submit large number of jobs using a loop or if you want to run a series of jobs that runs after completion of set of jobs using the same command. You can also schedule the job to start running on a predefined time as well. In this tutorial, we will explain how to submit jobs that runs depending on status of the previously submitted jobs or schedule a bunch of jobs to run one after the other.

To understand the dependency feature take a look at the `-d, --dependency` section in the `man` page of the `sbatch` command


Once you submit a job, using that job ID, you can submit dependency jobs. For eg.
<pre>
sbatch first_job.slurm
</pre>
You will get the job id
<pre>
854.condo
</pre>
Next, you can submit a job that only runs after successful completion of the first job as follows:
<pre>
sbatch --dependency=afterok:854 second_job.slurm
</pre>

The format here is 
<pre>
sbatch --dependency=type:job_id jobfile
</pre>

If the job requires more than one job to be completed before it is executed, you can supply all the jobids using `,` separator
<pre>
sbatch --dependency=type:job_id,job_id,job_id jobfile
</pre>

You can also set the job to run if any one of the job ids compltes successfully using a `?` separator
<pre>
sbatch --dependency=type:job_id?job_id?job_id jobfile
</pre>

The other dependencies that can be used for<blockcode><type:job_id></blockcode> are as follows:

<table>
<thead><tr><th>Argument</th><th>Description</th></tr></thead><tbody>
 <tr><td>after</td><td>This job can begin execution after the specified jobs have begun execution</td></tr>
 <tr><td>afterany</td><td>This job can begin execution after the specified jobs have terminated.</td></tr>
 <tr><td>aftercorr</td><td>A task of this job array can begin execution after the corresponding task ID in the specified job has completed successfully</td></tr>
 <tr><td>afternotok</td><td>This job can begin execution after the specified jobs have terminated in some failed state</td></tr>
 <tr><td>afterok</td><td>This job can begin execution after the specified jobs have successfully executed</td></tr>
 <tr><td>singleton</td><td>This job can begin execution after any previously launched jobs sharing the same job name and user have terminated</td></tr>
</tbody></table>
 

