#!/bin/bash

set -e
set -x  

git config --global user.name "Release Bot"
git config --global user.email "tabassum_kakar@hms.harvard.edu"

current_branch=$(git rev-parse --abbrev-ref HEAD)

echo "Fetching tags and ensuring full Git history..."
git fetch --unshallow 2>/dev/null || true  # Handle non-shallow repos gracefully
git fetch --tags

latest_tag=$(git tag | sort -V | tail -n 1)

echo "Latest tag found: $latest_tag"

# Function to get the next version tag
get_next_version_tag() {
  version_number="${latest_tag#v}"  # Remove the "v" prefix
  IFS='.' read -r major minor patch <<< "$version_number"
  
  new_patch=$((patch + 1))
  new_tag="v$major.$minor.$new_patch"

  echo "$new_tag"
}

# Function to generate release notes
generate_release_notes() {
  release_notes="Docker Version Changes in this release:\n"
  for change in "${version_changes[@]}"; do
    release_notes+="$change\n"
  done
  echo -e "$release_notes"
}

# Get the most recent tag
latest_tag=$(git tag | sort -V | tail -n 1)

# if [[ -z "$latest_tag" ]]; then
#   echo "No tags found. Starting with v0.0.1."
#   latest_tag="v0.0.1"
# fi

echo "Latest tag found: $latest_tag"

# Ensure HEAD is valid
if ! git rev-parse HEAD >/dev/null 2>&1; then
  echo "Error: HEAD is not pointing to a valid commit."
  exit 1
fi

# Get the list of commits since the last release
commits_since_last_tag=$(git log "$latest_tag"..HEAD --reverse --pretty=format:"%H")

for commit_hash in $commits_since_last_tag; do
  git checkout "$commit_hash" > /dev/null 2>&1

  echo "Checking commit $commit_hash..."

  version_changes=()
  version_changed=false

  # Iterate over all VERSION files in the containers directory
  while IFS= read -r version_file; do
    current_version=$(cat "$version_file")
    container_name="${version_file//containers\//}"
    container_name="${container_name%/VERSION}"

    echo "Version found in $version_file: $current_version"

    # For the first commit, there's no previous commit to compare to.
    if git cat-file commit "$commit_hash" >/dev/null 2>&1 && git show "$commit_hash^:$version_file" > /dev/null 2>&1; then
      previous_version=$(git show "$commit_hash^:$version_file" 2>/dev/null)
      
      if [[ "$previous_version" != "$current_version" ]]; then
        echo "Version change detected in $version_file:"
        echo "Previous version: $previous_version"
        echo "Current version:  $current_version"

        version_changed=true
        if [[ "$previous_version" == "" ]]; then
          version_changes+=("$container_name: Created v0.0.1")
        else
          version_changes+=("$container_name: $previous_version -> $current_version")
        fi
      fi
    else
      echo "First commit or no previous version found for $version_file."
      version_changed=true
      version_changes+=("$container_name: Created v0.0.1")
    fi
  done < <(find containers/* -type f -name "VERSION")

  # If any version change was detected, create a tag and release
  if $version_changed; then
    next_tag=$(get_next_version_tag)

    commit_message=$(git log -1 --format=%B "$commit_hash")

  if git tag -l | grep -q "$next_tag"; then
      echo "Tag $next_tag exists locally, deleting..."
      git tag -d "$next_tag"
    else
      echo "Tag $next_tag does not exist locally, skipping delete."
    fi

    # Create a tag
    echo "Creating tag: $next_tag"
    git tag -a "$next_tag" -m "Version changes in this release: ${version_changes[*]}" 

    # Generate release notes
    release_notes=$(generate_release_notes)

    # Create GitHub release with the generated release notes using GitHub
    echo "Creating GitHub release for tag: $next_tag"

    # Make the API request to create the release
    curl -X POST \
      -H "Authorization: token $RELEASE_TOKEN" \
      -d @- \
      https://api.github.com/repos/hubmapconsortium/portal-containers/releases <<EOF
{
  "tag_name": "$next_tag",
  "name": "Release $next_tag",
  "body": "$release_notes",
  "draft": false,
  "prerelease": false
}
EOF

  else
    echo "No version changes detected in this commit. Skipping tag creation."
  fi
done

echo "Returning to branch $current_branch..."
git checkout "$current_branch" > /dev/null 2>&1

echo "Version tagging completed."
