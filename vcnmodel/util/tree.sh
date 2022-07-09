#!/bin/bash

SEDMAGIC='s;[^/]*/;|-- ;g;s;-- |;   |;g'

if [ "$#" -gt 0 ] ; then
   dirlist="$@"
else
   dirlist="."
fi

for x in $dirlist; do
     find "$x" -print | sed -e "$SEDMAGIC"
done
