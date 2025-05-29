# GitHub Actions Workflows

## Active Workflows

- `deploy.yml`: Deploys the Jekyll website to GitHub Pages using Ruby/Jekyll.

## Disabled Workflows

The following workflows have been disabled to prevent CI failure due to Docker authentication issues:

- `deploy-image.yml.disabled`: This workflow was attempting to build and push a Docker image to Docker Hub.
- `deploy-docker-tag.yml.disabled`: This workflow was attempting to tag and push Docker images when new version tags were created.

These workflows were disabled because they required Docker Hub authentication credentials (`DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN`) that weren't configured in the repository secrets.

If you need to use these Docker-based workflows in the future:

1. Create Docker Hub credentials
2. Add them as repository secrets in GitHub:
   - Go to your repository → Settings → Secrets and variables → Actions
   - Add `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN` secrets
3. Rename the files back to remove the `.disabled` extension
