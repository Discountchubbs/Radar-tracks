#!/bin/bash

# Abort any existing cherry-pick
git cherry-pick --abort

# Start the cherry-pick
git cherry-pick -n 9ffbbf8022be64213af2b83ff6208b4687573671

# Force remove all conflicting files
git ls-files -u | cut -f 2 | sort -u | while read file; do
    git rm -f "$file"
    git checkout CHERRY_HEAD -- "$file" 2>/dev/null || true
done

# Add all changes
git add -A

# Complete the cherry-pick
git commit -C 9ffbbf8022be64213af2b83ff6208b4687573671
EOF
