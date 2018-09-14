#! /bin/bash

if [[ $# != 1 ]]; then
  echo "Usage: ${0} prefix"
  echo "  prefix: 文件名前缀，包括输入和输出"
  echo "将前面的程序产生的一系列网格 Open DX 的文件，产生 prefix.dx 文件将这些网格文件组合到一起用于显示。"
  exit 1
fi

PREFIX=$1
rm -f ${PREFIX}.dx
tmpfile=`mktemp ${PREFIX}.XXXXXX`
echo "object \"group\" class group" >${tmpfile}
value=0
for dxi in `ls ${PREFIX}?*.dx`; do 
  echo -n "${dxi} "
  echo -n "  member \"`basename ${dxi}`\"" >>${tmpfile}
  echo "  value file \"${dxi}\"" >>${tmpfile}
  pnt=`grep "^object 1" ${dxi} | cut -d ' ' -f 12`
  dxfile=`mktemp ${dxi}.XXXXXX`
  sed -e "s/\(attribute \"ref\" string \"positions\"\)/\1\n\nobject 3 class constantarray type float rank 1 shape 1 item ${pnt} data follows\n${value} \nattribute \"dep\" string \"positions\"/g" ${dxi} | sed -e "s/\(component \"connections\" value 2\)/\1\ncomponent \"data\" value 3/g" >${dxfile}
  mv ${dxfile} ${dxi}
  let "value += 1"
done
echo "done!"
mv ${tmpfile} ${PREFIX}.dx

exit 0

