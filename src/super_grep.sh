#! /bin/bash 


if [ $# -eq 2 ];then
  grep "$1" . -inr --include *.$2
elif [ $# -eq 3 ];then
  grep "$1" . -inr --include *.{$2,$3}
elif [ $# -eq 4 ];then
  grep "$1" . -inr --include *.{$2,$3,$4}
elif [ $# -eq 5 ];then
  grep "$1" . -inr --include *.{$2,$3,$4,$5}
else
  echo "usage : ./super_grep $1 $n"
  exit
fi


