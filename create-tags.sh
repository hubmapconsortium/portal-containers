#!/bin/bash

# Optional dry-run flag
DRY_RUN=false

# Capture the current branch name to return back to it later
current_branch=$(git rev-parse --abbrev-ref HEAD)

# Check if --dry-run flag is provided
if [[ "$1" == "--dry-run" ]]; then
  DRY_RUN=true
fi

# Function to get the next version tag, incrementing the patch version
get_next_version_tag() {
  # Get the most recent tag (sorted in version order)
  latest_tag=$(git tag | sort -V | tail -n 1)

  # If there are no tags, start from v0.0.1
  if [[ -z "$latest_tag" ]]; then
    echo "v0.0.1"
    return
  fi

  # Extract the version number and increment it
  version_number="${latest_tag#v}"  # Remove the "v" prefix
  IFS='.' read -r major minor patch <<< "$version_number"
  
  # Increment the patch version (you can change this to major or minor depending on your needs)
  new_patch=$((patch + 1))
  
  # Construct the new tag name
  new_tag="v$major.$minor.$new_patch"

  echo "$new_tag"
}

# Iterate over all commits in reverse order (from oldest to newest)
git log --reverse --pretty=format:"%H" | while read commit_hash; do
  # Check out the commit (to inspect the state at that commit)
  git checkout "$commit_hash" > /dev/null 2>&1

  echo "Checking commit $commit_hash..."

  # Initialize variables to track changes
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
    # Get the next tag (incremented version)
    next_tag=$(get_next_version_tag)

    # Get the commit message for this commit
    commit_message=$(git log -1 --format=%B "$commit_hash")

    if [ "$DRY_RUN" = false ]; then
      # Create the tag
      echo "Creating tag: $next_tag"
      git tag -a "$next_tag" -m "Version $next_tag: $commit_message"

      # Add git notes with version changes
      git notes add -m "Version changes in this commit: ${version_changes[*]}"
    else
      echo "[DRY-RUN] Would create tag: $next_tag"
      echo "[DRY-RUN] Would add git note: Version changes in this commit: ${version_changes[*]}"
    fi
  else
    echo "No version changes detected in this commit. Skipping tag creation."
  fi
done

# Checkout back to the latest branch
echo "Returning to branch $current_branch..."
git checkout "$current_branch" > /dev/null 2>&1

echo "Version tagging completed."
