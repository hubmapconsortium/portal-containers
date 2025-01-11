#!/bin/bash

# Log file where all output messages will be written
LOG_FILE="version_tagging.log"
> "$LOG_FILE"  # Clear the log file at the start of each run

echo "Starting version tagging process..." | tee -a "$LOG_FILE"

# Capture the current branch name to return back to it later
current_branch=$(git rev-parse --abbrev-ref HEAD)

# Stash any uncommitted changes so we don't lose them
git stash push -m "Temporary stash before version tagging" > /dev/null 2>&1

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
  
  # Increment the patch version
  new_patch=$((patch + 1))
  
  # Construct the new tag name
  new_tag="v$major.$minor.$new_patch"

  echo "$new_tag"
}

# Store versions from the previous commit
declare -A previous_versions

# Iterate over the first 20 commits in reverse order (from oldest to newest)
git log --reverse --pretty=format:"%H" | head -n 20 | while read commit_hash; do
  # Check out the commit (to inspect the state at that commit)
  git checkout "$commit_hash" > /dev/null 2>&1

  echo "Processing commit: $commit_hash" | tee -a "$LOG_FILE"

  # Flag to check if any version change is detected
  version_changed=false
  updated_containers=""

  # Check if the containers directory exists in this commit
  if [ ! -d "containers" ]; then
    echo "No containers directory found in this commit. Skipping version check." | tee -a "$LOG_FILE"
    continue
  fi

  # Check if there are any VERSION files in the containers directory
  version_files=$(find containers/* -type f -name "VERSION")
  
  # If there are no VERSION files, skip this commit
  if [[ -z "$version_files" ]]; then
    echo "No VERSION files found in this commit. Skipping tag creation." | tee -a "$LOG_FILE"
    continue
  fi

  # Iterate over all VERSION files in the containers directory
  find containers/* -type f -name "VERSION" | while read version_file; do
    container_name=$(basename $(dirname "$version_file"))
    current_version=$(cat "$version_file")

    echo "Version found in $version_file for container $container_name: $current_version" | tee -a "$LOG_FILE"

    # If the version has changed since the previous commit, create a tag
    # Special case: Treat version 0.0.1 as a version change if it's the first time seeing it
    if [[ "$current_version" == "0.0.1" && -z "${previous_versions[$container_name]}" ]]; then
      echo "First version found for $container_name: $current_version. Treating as a version change." | tee -a "$LOG_FILE"
      version_changed=true
      updated_containers="$updated_containers $container_name"
    elif [[ -n "${previous_versions[$container_name]}" ]] && [[ "${previous_versions[$container_name]}" != "$current_version" ]]; then
      echo "Version change detected for $container_name: ${previous_versions[$container_name]} -> $current_version" | tee -a "$LOG_FILE"
      version_changed=true
      updated_containers="$updated_containers $container_name"
    fi

    # Store the current version for the next commit comparison
    previous_versions[$container_name]="$current_version"
  done

  # If a version change was detected, create a new tag for the commit
  if $version_changed; then
    next_tag=$(get_next_version_tag)
    commit_message=$(git log -1 --format=%B "$commit_hash")

    # Create the tag with commit message as release notes
    echo "Creating tag: $next_tag for commit $commit_hash, containers updated: $updated_containers" | tee -a "$LOG_FILE"
    git tag -a "$next_tag" -m "Version $next_tag for containers: $updated_containers. $commit_message"

  else
    echo "No version changes detected in this commit. Skipping tag creation." | tee -a "$LOG_FILE"
  fi

done

# Checkout back to the branch in a clean state without changing anything in the working directory
git checkout "$current_branch" > /dev/null 2>&1

# Apply the stashed changes back to the working directory
git stash pop > /dev/null 2>&1

# Return to the branch in a clean state
echo "Returning to branch $current_branch with uncommitted changes restored..." | tee -a "$LOG_FILE"

echo "Version tagging completed." | tee -a "$LOG_FILE"
