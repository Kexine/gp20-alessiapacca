#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR/../..

IFS=$'\n'
FILES=($(pkgutil --files com.intel.embree-2.17.4 | grep -e 'opt/local/include/\|opt/local/lib/\|Applications/Embree2' | tail -r))
unset IFS

# exit if no files found
if [ ${#FILES[@]} -eq 0 ]; then
  printf "Embree 2.17.4 not installed!\n"
  exit
fi

# first print all files that would get removed
echo Uninstalling Embree 2.17.4 will remove the following files:
PWD=`pwd`
if [ "$PWD" != "/" ]; then
    PWD=$PWD/
fi
for f in "${FILES[@]}"; do
    printf "  %s%s\n" $PWD "$f"
done

echo "Do you wish to uninstall Embree 2.17.4 by removing these files?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done

# now remove files
echo Uninstalling Embree 2.17.4 ...
for f in "${FILES[@]}"; do
    sudo /bin/rm -vd "$f"
done

sudo /usr/sbin/pkgutil --forget com.intel.embree-2.17.4
