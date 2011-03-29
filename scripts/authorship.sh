#!/bin/bash

git log --numstat | gawk '/^Author:/ {
   author           = $2
   commits[author] += 1
   commits["tot"]  += 1
}

/^[0-9]+ +[0-9]+ +vendor|framework/ { next }

/^[0-9]/ {
   more[author] += $1
   less[author] += $2
   file[author] += 1

   more["tot"]  += $1
   less["tot"]  += $2
   file["tot"]  += 1
}

END {
   for (author in commits) {
      if (author != "tot") {
         more[author]    = more[author] / more["tot"] * 100
         less[author]    = less[author] / less["tot"] * 100
         file[author]    = file[author] / file["tot"] * 100
         commits[author] = commits[author] / commits["tot"] * 100

         printf "%s:\n  insertions: %.0f%%\n  deletions: %.0f%%\n  files: %.0f%%\n  commits: %.0f%%\n", author, more[author], less[author], file[author], commits[author]
      }
   }
}'
