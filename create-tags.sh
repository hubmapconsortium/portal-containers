#!/bin/bash

current_branch=$(git rev-parse --abbrev-ref HEAD)

get_next_version_tag() {

  latest_tag=$(git tag | sort -V | tail -n 1)


  if [[ -z "$latest_tag" ]]; then
    echo "v0.0.1"
    return
  fi

 
  version_number="${latest_tag#v}"  # Remove the "v" prefix
  IFS='.' read -r major minor patch <<< "$version_number"
  
  new_patch=$((patch + 1))
  
  new_tag="v$major.$minor.$new_patch"

  echo "$new_tag"
}


git log --reverse --pretty=format:"%H" | while read commit_hash; do
  git checkout "$commit_hash" > /dev/null 2>&1

  echo "Checking commit $commit_hash..."

  version_changes=()
  version_changed=false

  # Iterate over all VERSION files in the containers directory
  while IFS= read -r version_file; do
    current_version=$(cat "$version_file")
    echo "Version found in $version_file: $current_version"

    # For the first commit, there's no previous commit to compare to.
    if git cat-file commit "$commit_hash" >/dev/null 2>&1 && git show "$commit_hash^:$version_file" > /dev/null 2>&1; then
      previous_version=$(git show "$commit_hash^:$version_file" 2>/dev/null)
      
      if [[ "$previous_version" != "$current_version" ]]; then
        echo "Version change detected in $version_file:"
        echo "Previous version: $previous_version"
        echo "Current version:  $current_version"

        version_changed=true
        version_changes+=("$version_file: $previous_version -> $current_version")
      fi
    else
      echo "First commit or no previous version found for $version_file."
      version_changed=true
      version_changes+=("$version_file: (initial) -> $current_version")
    fi
  done < <(find containers/* -type f -name "VERSION")

  # If any version change was detected, create a tag and add git notes
  if $version_changed; then
    next_tag=$(get_next_version_tag)

    commit_message=$(git log -1 --format=%B "$commit_hash")

    # if [ "$DRY_RUN" = false ]; then
    echo "Creating tag: $next_tag"
    git tag -a "$next_tag" -m "Version $next_tag: $commit_message"

    git notes add -f -m "Version changes in this commit: ${version_changes[*]}"
    # else
    #   echo "[DRY-RUN] Would create tag: $next_tag"
    #   echo "[DRY-RUN] Would add git note: Version changes in this commit: ${version_changes[*]}"
    # fi
  else
    echo "No version changes detected in this commit. Skipping tag creation."
  fi
done

# Checkout back to the latest branch
echo "Returning to branch $current_branch..."
git checkout "$current_branch" > /dev/null 2>&1

echo "Version tagging completed."
