USER=csw
HOST=mauritius.bchem.washington.edu

all:
	(cd ..; zip -r --verbose profiling.zip profiling shared)
	scp ../profiling.zip ${USER}@${HOST}:ShinyApps/
	ssh -Y ${USER}@${HOST} '(cd ShinyApps; unzip profiling.zip)'
