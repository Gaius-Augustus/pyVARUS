#!/bin/bash
# Build and push the pyVARUS image from the checked-out source tree.
# Run on a machine with Docker and `docker login` done (e.g. greif14).
#
#   docker/build_push.sh                 # tags v<pyproject version> and latest
#   TAG=v2.0.0 docker/build_push.sh      # explicit version tag (plus latest)
#   PUSH=0 docker/build_push.sh          # build only
#   DOCKER="sudo docker" docker/build_push.sh
#
# Refuses to build from a dirty working tree so the image always matches a
# commit; set ALLOW_DIRTY=1 to override for a local test build.
#
# Users are told to pull :latest (README). The version tag is kept so
# pipelines such as BRAKER4 can pin a reproducible image.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

DOCKER=${DOCKER:-docker}
IMAGE=${IMAGE:-gaiusaugustus/pyvarus}
VERSION=$(sed -n 's/^version = "\(.*\)"/\1/p' pyproject.toml | head -1)
TAG=${TAG:-v$VERSION}
PUSH=${PUSH:-1}
COMMIT=$(git rev-parse HEAD)

if [ -n "$(git status --porcelain -- varus pyproject.toml README.md LICENSE docker)" ] \
        && [ "${ALLOW_DIRTY:-0}" != 1 ]; then
    echo "Working tree has uncommitted changes in files that go into the image." >&2
    echo "Commit them first, or set ALLOW_DIRTY=1 for a throwaway test build." >&2
    exit 1
fi

echo "[$(date)] Building $IMAGE:$TAG from $(git rev-parse --abbrev-ref HEAD)@${COMMIT:0:7}"
$DOCKER build -f docker/Dockerfile \
    --build-arg VARUS_GIT_COMMIT="$COMMIT" \
    -t "$IMAGE:$TAG" -t "$IMAGE:latest" .

if [ "$PUSH" = 1 ]; then
    echo "[$(date)] Pushing $IMAGE:$TAG"
    $DOCKER push "$IMAGE:$TAG"
    $DOCKER push "$IMAGE:latest"
fi

echo "[$(date)] Done."
echo ""
echo "Pull with:"
echo "  singularity pull pyvarus.sif docker://$IMAGE:latest"
