ORG=$1
NEW=$2
WHERE=$3

find $WHERE -type f -exec sed -i "s/${ORG}/${NEW}/g" {} \;

find $WHERE -type f -exec rename $ORG $NEW {} \;