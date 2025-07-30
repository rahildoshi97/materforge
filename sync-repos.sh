#!/bin/bash
echo "Fetching from both remotes..."
git fetch github
git fetch origin

echo "Syncing GitHub main -> GitLab master..."
git push origin github/main:master

echo "Syncing GitLab master -> GitHub main..."
git push github origin/master:main

echo "Sync complete!"
