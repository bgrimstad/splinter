#!/bin/bash

ROOT=$(pwd)
OS="unknown"
COMPILER="unknown"
NPROC="unknown"

# Try to find CMakeLists.txt in this directory or its parent
SPLINTER_DIR="./"
if [ ! -f $SPLINTER_DIR/CMakeLists.txt ]; then
       SPLINTER_DIR="../"
       if [ ! -f $SPLINTER_DIR/CMakeLists.txt ]; then
               echo "Error: Unable to locate CMakeLists.txt!"
               exit 1
       fi
fi

# Defaults
CMAKE_CMD="cmake"

#MinGW config
MINGW_32_BIT="/C/mingw-w64/i686-4.9.2-posix-dwarf-rt_v4-rev3/mingw32/bin"
MINGW_64_BIT="/C/mingw-w64/x86_64-4.9.2-posix-seh-rt_v4-rev3/mingw64/bin"

# MSVC config
#MSBUILD_DIR="/C/Program Files (x86)/MSBuild/12.0/Bin"
#VCVARSALL_DIR="/C/Program Files (x86)/Microsoft Visual Studio 12.0/VC"
#MSVC_GENERATOR="Visual Studio 12 2013"
MSBUILD_DIR="/C/Program Files (x86)/MSBuild/14.0/Bin"
VCVARSALL_DIR="/C/Program Files (x86)/Microsoft Visual Studio 14.0/VC"
MSVC_GENERATOR="Visual Studio 14 2015"

# Capture the command argument for use in help messages
COMMAND=$0

# Thanks to http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash#14203146
# for this great command line argument parsing algorithm
# Use > 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
while [[ $# > 0 ]]
do
key="$1"

case $key in
	-m|--mingw-binary-dir)
	PATH="$2:$PATH"
	shift # past argument
	;;
	-c|--cmake-binary-dir)
	PATH="$2:$PATH"
	shift # past argument
	;;
	-vc|--vcvarsall-dir)
	VCVARSALL_DIR="$2"
	shift # past argument
	;;
	-mb|--msbuild-dir)
	MSBUILD_DIR="$2"
	shift # past argument
	;;
	*)
	# No preceding: path to SPLINTER
	SPLINTER_DIR=$1
	;;
esac
shift # past argument or value
done

# Make sure SPLINTER_DIR is an absolute path
if [[ $SPLINTER_DIR != /* ]]; then
	SPLINTER_DIR="$(pwd)/$SPLINTER_DIR"
fi

BUILD_ROOT=$SPLINTER_DIR/build

# Read the SPLINTER version from version.txt
SPLINTER_VERSION=$(cat $SPLINTER_DIR/version.txt)

# Check that we can find CMake
if [[ $(which cmake) == "" ]]; then
	echo "Error: Can't find CMake, make sure it is in your PATH environment variable"
	echo "and try again!"
	echo "If you don't want to add CMake to your PATH, you can specify the path to it with:"
	echo "$COMMAND -c /path/to/cmake/binary/directory"
	echo "Note that on Windows, the path \"C:/Program Files (x86)/CMake/bin\" has to be written as \"/c/Program Files (x86)/CMake/bin\""
	exit 1
fi

function get_commit_id {
	COMMIT_ID=$(git -C "$SPLINTER_DIR" log -n 1 --pretty=format:"%H")
}

function update_commit_id {
	get_commit_id
	echo $COMMIT_ID > "$BUILD_ROOT/$OS/$COMPILER/commit_id"
}

function update_compiler_version {
	echo $COMPILER_VERSION > "$BUILD_ROOT/$OS/$COMPILER/compiler_version"
}

function copy_header_files {	
	cp -r $SPLINTER_DIR/include $BUILD_ROOT/$OS/$COMPILER/
	cp -r $SPLINTER_DIR/thirdparty/Eigen $BUILD_ROOT/$OS/$COMPILER/include
}

function build_gcc_clang {
	ARCH=$1
	COMPILER=$2

	mkdir -p $BUILD_ROOT/$OS/$COMPILER/$ARCH
	mkdir -p $BUILD_ROOT/.build/$COMPILER/$ARCH
	cd $BUILD_ROOT/.build/$COMPILER/$ARCH
	
	rm CMakeCache.txt
	echo "Building SPLINTER for $ARCH with $COMPILER"
	"$CMAKE_CMD" "$SPLINTER_DIR" -DCMAKE_BUILD_TYPE=release -DARCH=$ARCH -G "Unix Makefiles" -DCMAKE_MAKE_PROGRAM="$MAKE_CMD"
	"$MAKE_CMD" -j$NPROC
}

function build_linux {
	echo "Building for Linux"
	OS=linux
	
	MAKE_CMD=$(which make)
	NPROC=$(nproc)
	
	GPP=$(which g++)
	if [[ $GPP != "" ]]; then
		export CXX=$GPP
		COMPILER=gcc
		COMPILER_VERSION=$($CXX -dumpversion)
		
		build_gcc_clang x86 $COMPILER
		cp libsplinter-$SPLINTER_VERSION.so libsplinter-static-$SPLINTER_VERSION.a "$BUILD_ROOT/$OS/$COMPILER/$ARCH"
		# MatLab for Linux only exists as 64bit, so we don't need this
#		"$MAKE_CMD" install
#		cp -r splinter-matlab $BUILD_ROOT
		
		build_gcc_clang x86-64 $COMPILER
		cp libsplinter-$SPLINTER_VERSION.so libsplinter-static-$SPLINTER_VERSION.a "$BUILD_ROOT/$OS/$COMPILER/$ARCH"
		"$MAKE_CMD" install
		cp -r splinter-matlab $BUILD_ROOT
		
		copy_header_files
		
		# Write down the commit id this was compiled from
		update_commit_id
		update_compiler_version
	fi
	
	CLANG=$(which clang++-3.5)
	if [[ $CLANG != "" ]]; then
		export CXX=$CLANG
		COMPILER=clang
		COMPILER_VERSION=$($CXX -dumpversion)
		
		build_gcc_clang x86 $COMPILER
		cp libsplinter-$SPLINTER_VERSION.so libsplinter-static-$SPLINTER_VERSION.a "$BUILD_ROOT/$OS/$COMPILER/$ARCH"
		
		build_gcc_clang x86-64 $COMPILER
		cp libsplinter-$SPLINTER_VERSION.so libsplinter-static-$SPLINTER_VERSION.a "$BUILD_ROOT/$OS/$COMPILER/$ARCH"
		
		copy_header_files
		
		# Write down the commit id this was compiled from
		update_commit_id
		update_compiler_version
	fi
}

function build_msvc {
	ARCH=$1
	COMPILER=$2
	
	mkdir -p $BUILD_ROOT/$OS/$COMPILER/$ARCH
	mkdir -p $BUILD_ROOT/.build/$COMPILER/$ARCH
	cd $BUILD_ROOT/.build/$COMPILER/$ARCH
	
	# Need this so msbuild.exe can find the project file
#	export PATH="$ROOT/build/$COMPILER/$ARCH/:$PATH"
	rm CMakeCache.txt
	
	if [[ $ARCH == "x86" ]]; then
		cmd "/C vcvarsall.bat x86"
		GENERATOR=$MSVC_GENERATOR
	elif [[ $ARCH == "x86-64" ]]; then
		cmd "/C vcvarsall.bat x64"
		GENERATOR="$MSVC_GENERATOR Win64"
	else
		echo "Error: Unknown architecture given to build_msvc: $ARCH"
		exit 1
	fi
	
	"$CMAKE_CMD" "$SPLINTER_DIR" -DCMAKE_BUILD_TYPE=Release -DARCH=$ARCH -G "$GENERATOR"
	
	"$MSBUILD" ALL_BUILD.vcxproj -p:Configuration=Release -maxcpucount:$NPROC
	
	# Install
	mkdir -p "$BUILD_ROOT/splinter-matlab/lib/$OS/$ARCH/"
	mkdir -p "$BUILD_ROOT/$OS/$COMPILER/$ARCH"
	
	cp "Release/splinter-matlab-$SPLINTER_VERSION.dll" "$BUILD_ROOT/splinter-matlab/lib/$OS/$ARCH/"
	cp "Release/splinter-$SPLINTER_VERSION.dll" "$BUILD_ROOT/$OS/$COMPILER/$ARCH"
	cp "Release/splinter-static-$SPLINTER_VERSION.lib" "$BUILD_ROOT/$OS/$COMPILER/$ARCH"

	copy_header_files
}

function build_windows {
	echo "Building for Windows"
	OS=windows
	
	export PATH="$MSBUILD_DIR:$PATH"
	export PATH="$VCVARSALL_DIR:$PATH"
	
	# Get number of processors for use with -maxcpucount
	NPROC_STRING=$(cmd "/C echo %NUMBER_OF_PROCESSORS%")
	NPROC="${NPROC_STRING//[!0-9]/}"
	
	# First build with MinGW if it is installed and in PATH
	GPP="g++"
	MAKE_CMD="mingw32-make"
	if [[ $GPP != "" && $MAKE_CMD != "" ]]; then
		export CXX=$GPP
		COMPILER=gcc
		COMPILER_VERSION=$($CXX -dumpversion)

		if [[ $MINGW_64_BIT != "" ]]; then
			export PATH="$MINGW_32_BIT:$PATH"
			build_gcc_clang x86 $COMPILER
		fi
		
		cp libsplinter-$SPLINTER_VERSION.dll libsplinter-static-$SPLINTER_VERSION.a "$BUILD_ROOT/$OS/$COMPILER/$ARCH"

		if [[ $MINGW_64_BIT != "" ]]; then
			export PATH="$MINGW_64_BIT:$PATH"
			build_gcc_clang x86-64 $COMPILER
		fi
		
		copy_header_files
		
		# Write down the commit id this was compiled from
		update_commit_id
		update_compiler_version
	fi

	MSBUILD=$(which msbuild.exe)
	if [[ $MSBUILD != "" ]]; then
		COMPILER=msvc
		COMPILER_VERSION=$(msbuild.exe "-version" | grep '^[[:digit:]]\+.[[:digit:]]\+.[[:digit:]]\+.[[:digit:]]\+$')
		
		last_compiled_commit_id=$(cat $BUILD_ROOT/$OS/$COMPILER/commit_id)
		get_commit_id
		if [[ $last_compiled_commit_id == $COMMIT_ID ]]; then
			echo "No new commits since last compile with $COMPILER, skipping."
			
		else
			build_msvc "x86" $COMPILER
			build_msvc "x86-64" $COMPILER
			
			# Write down the commit id this was compiled from
			update_commit_id
			update_compiler_version
		fi
	fi
}


mkdir -p $BUILD_ROOT # -p to avoid error message when it already exists
cd $BUILD_ROOT

PLATFORM=$(uname)
if [[ $PLATFORM == MINGW* ]]; then
	build_windows
	
elif [[ $PLATFORM == Linux ]]; then
	build_linux
	
else
	echo "Unknown platform: $PLATFORM"
fi

cd $BUILD_ROOT
# Check that all commit ids are the same
# If they are we can make a release
# TODO: Add osx
OSES="windows
linux"
COMMIT_ID=""
for os_dir in $OSES
do
	if [[ ! -d $os_dir ]]; then
		echo "Cannot make release because an OS directory ($os_dir) is missing."
		exit 1
	fi
	
	for compiler in $(ls $BUILD_ROOT/.build/$os_dir)
	do
		if [[ $COMMIT_ID == "" ]]; then
			COMMIT_ID=$(cat $BUILD_ROOT/.build/$os_dir/$compiler/commit_id)
		else
			if [[ $(cat $BUILD_ROOT/.build/$os_dir/$compiler/commit_id) != $COMMIT_ID ]]; then
				echo "Commit id mismatch, $os_dir/$compiler differs from previous."
				echo "Cannot make release."
				exit 1
			fi
		fi
	done
done

echo "All builds were built from the same commit, proceeding to make release."

# If tar is installed, and all commit ids are the same,
# then we make a release
TAR=$(which tar)
ZIP=$(which zip)
if [[ $TAR == ""  || $ZIP == "" ]]; then
	echo "Error: Missing either tar or zip, need both to create release."
	exit 1
fi

mkdir -p $BUILD_ROOT/releases
rm $BUILD_ROOT/releases/*
for os_dir in $OSES
do
	cd $BUILD_ROOT/.build/$os_dir
	for compiler_dir in $(echo */) # echo */ gives us a list of the directories
	do
		compiler_name=${compiler_dir%?} # compiler_dir includes the last /, remove it.
		cd $BUILD_ROOT/.build/$os_dir/$compiler_dir
		files=""
		for arch in $(echo */)
		do
			files="$arch $files"
		done

		filename=$os_dir"_"$compiler_name$(cat compiler_version)
		full_filename=$BUILD_ROOT/releases/$filename

		OLDWD=$(pwd)
		cd $BUILD_ROOT/.build/$os_dir/$compiler_dir

		echo "Creating archive $filename.tar.gz"
		$TAR -czf $full_filename.tar.gz $files > /dev/null

		echo "Creating archive $filename.zip"
		$ZIP -r $full_filename $files > /dev/null
		cd $OLDWD
	done
done

# Make an archive of splinter-matlab
filename="splinter-matlab"
full_filename="$BUILD_ROOT/releases/$filename"
files="splinter-matlab"
cd $BUILD_ROOT
echo "Creating archive $filename.tar.gz"
$TAR -czf $full_filename.tar.gz $files > /dev/null

echo "Creating archive $filename.zip"
$ZIP -r $full_filename $files > /dev/null
