# Let's run the scheduler and science on a few runs
#tail -n+4 generate.sh | xargs -n 3 -I'{}' -P3 bash -c '{}'

'rm -r minion_1016/sched'
'rm -r minion_1012/sched'
'rm -r kraken_1042/sched'
'rm -r minion_1016/sci'
'rm -r minion_1012/sci'
'rm -r kraken_1042/sci'

'schedulerValidation.py  ~/Scratch/Opsim_sqlites/minion_1016_sqlite.db   --outDir minion_1016/sched'
'schedulerValidation.py  ~/Scratch/Opsim_sqlites/minion_1012_sqlite.db	--outDir minion_1012/sched'
'schedulerValidation.py  ~/Scratch/Opsim_sqlites/kraken_1042_sqlite.db	--outDir kraken_1042/sched'
'sciencePerformance.py   ~/Scratch/Opsim_sqlites/minion_1016_sqlite.db	--outDir minion_1016/sci  '
'sciencePerformance.py   ~/Scratch/Opsim_sqlites/minion_1012_sqlite.db	--outDir minion_1012/sci  '
'sciencePerformance.py   ~/Scratch/Opsim_sqlites/kraken_1042_sqlite.db	--outDir kraken_1042/sci  '
