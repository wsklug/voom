#!/bin/csh -f
#  hex50L.cmd
#
#  UGE job for indent built Mon Aug  6 21:36:55 PDT 2012
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/campus/klug/Software/voom/src/Applications/graphene/hex50L.joblog.$JOB_ID
#$ -o /u/home/campus/klug/Software/voom/src/Applications/graphene/hex50L.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/campus/klug/Software/voom/src/Applications/graphene/indent
#  arguments       = hex50L -i 200.0 -m 1.0 -s 1.0 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#
#$ -l h_data=1024M,h_rt=24:00:00
#
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M klug@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#  Uncomment the next line to have your environment variables used by SGE
# #$ -V
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqidir    = /u/home/campus/klug/Software/voom/src/Applications/graphene
  set qqjob     = indent
  set qqodir    = /u/home/campus/klug/Software/voom/src/Applications/graphene
  cd     /u/home/campus/klug/Software/voom/src/Applications/graphene
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for indent built Mon Aug  6 21:36:55 PDT 2012"
  echo ""
  echo "  indent directory:"
  echo "    "/u/home/campus/klug/Software/voom/src/Applications/graphene
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "indent started on:   "` hostname -s `
  echo "indent started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load intel/11.1
  module load vtk
#
  echo indent "hex50L -i 200.0 -m 1.0 -s 1.0" \>\& hex50L.output.$JOB_ID
  echo ""
  time /u/home/campus/klug/Software/voom/src/Applications/graphene/indent hex50L -i 200.0 -m 1.0 -s 1.0 >& /u/home/campus/klug/Software/voom/src/Applications/graphene/hex50L.output.$JOB_ID
#
  echo ""
  echo "indent finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/campus/klug/Software/voom/src/Applications/graphene/hex50L.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
  if (`wc -l /u/home/campus/klug/Software/voom/src/Applications/graphene/hex50L.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
	head -50 /u/home/campus/klug/Software/voom/src/Applications/graphene/hex50L.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
	echo " "  >> /u/local/apps/queue.logs/job.log.serial
	tail -10 /u/home/campus/klug/Software/voom/src/Applications/graphene/hex50L.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
	cat /u/home/campus/klug/Software/voom/src/Applications/graphene/hex50L.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)
