#! /bin/bash

if [[ $# != 1 ]]; then
  echo "Usage: ${0} prefix"
  echo "  prefix: 文件名前缀，包括输入和输出"
  echo "将前面的程序产生的一系列有限元函数 Open DX 的文件，产生 prefix.dx 文件将这些文件组合到一起用于显示。"
  exit 1
fi

PREFIX=$1
tmpfile=`mktemp ${PREFIX}.XXXXXX`
echo "object \"group\" class group" >${tmpfile}
value=0
for dxi in `ls ${PREFIX}?*.dx`; do 
  echo -n "${dxi} "
  echo -n "  member \"`basename ${dxi}`\"" >>${tmpfile}
  echo "  value file \"${dxi}\"" >>${tmpfile}
  let "value += 1"
done
echo "done!"
mv ${tmpfile} ${PREFIX}.dx

exit 0

