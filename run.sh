# !/bin/bash

DDAR_ARGS=(
  --defs_file=$(pwd)/data/defs.txt \
  --rules_file=$(pwd)/data/rules.txt \
);

SEARCH_ARGS=(
  --batch_size=2
  --beam_size=2
  --search_depth=2
)

python src/alphageometry.py \
--alsologtostderr \
--problems_file=$(pwd)/examples/examples.txt \
--problem_name=orthocenter \
--mode=alphageometry \
"${DDAR_ARGS[@]}" \
"${SEARCH_ARGS[@]}"
