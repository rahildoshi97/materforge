#!/bin/bash
echo "Fetching from both remotes..."
git fetch github
git fetch origin

echo "Syncing GitHub main -> GitLab pwlf..."
git push origin github/main:pwlf

echo "Syncing GitLab master -> GitHub main..."
git push github origin/master:main

echo "Sync complete!"
