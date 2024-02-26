#!/bin/bash

# Get the list of all installed packages with their versions inside Docker container
#pip list --format=json > requirements.json

# Get the list of dependencies for pygmtsar
dependencies=$(pip3 show pygmtsar 2>/dev/null | grep Requires | cut -d: -f2 | tr -d ' ')

echo "RUN pip3 install \\"

# Loop through the dependencies to find their versions
for dep in ${dependencies//,/ }; do
  # Use jq to parse the JSON and find the version of the dependency
  version=$(jq --arg dep "$dep" -r '.[] | select(.name==$dep) | .version' requirements.json)
  if [[ ! -z "$version" ]]; then
    # Write the dependency and its version to requirements.txt
    echo "    $dep==$version \\"
  fi
done

echo "    pygmtsar"
