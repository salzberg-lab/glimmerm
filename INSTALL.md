# How to Compile the Code on Linux

```bash
# Go to the directory containing the Makefile
cd sources

# Compile the source code and create the executables
make

# Copy the executables to the bin/ directory
cp glimmerhmm ../bin/

# Go to the bin directory
cd ../bin/

# List the executable files with details
ls -l

# Add the absolute path of the bin directory to PATH for this shell session
export PATH="$(pwd -P):$PATH"
