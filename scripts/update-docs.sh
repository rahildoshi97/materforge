#!/bin/bash
set -e

echo "Building documentation locally..."
cd docs
rm -rf _build/
sphinx-build -b html . _build/html

echo "Documentation built successfully!"
echo "View locally at: file://$(pwd)/_build/html/index.html"
echo "Or run: cd _build/html && python -m http.server 8000"

echo "Checking for documentation changes..."
cd ..

# Add any new files and check if there are changes
git add docs/ update-docs.sh

if git diff --cached --quiet; then
    echo "No documentation changes to commit."
else
    echo "Committing documentation changes..."
    git commit -m "Update documentation"
    
    echo "Pushing to GitHub (this will trigger Read the Docs rebuild)..."
    git push github main
    
    echo "Syncing to GitLab..."
    ./sync-repos.sh
    
    echo "Done! Check https://materforge.readthedocs.io/ in a few minutes for updated online docs."
fi
