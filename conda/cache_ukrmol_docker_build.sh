#!/bin/bash
# ============================================================================ #
# UKRmol+ Docker Cache Extractor
# ============================================================================ #
# DESCRIPTION:
#   This script extracts build artifacts (downloaded source tarballs and 
#   compiled binaries) from an existing UKRmol+ Docker image and saves them 
#   to a local './ukrmol-cache' directory.
#
# PURPOSE:
#   To dramatically speed up future Docker builds. When you run 'docker build',
#   the Dockerfile detects the 'ukrmol-cache' folder and injects these files,
#   skipping the time-consuming download and compilation steps (e.g., GCC).
#
# LOGIC:
#   The script automatically detects the type of image provided:
#     1. Builder Stage (Preferred): Extracts downloads AND binaries from /build.
#     2. Final Stage (Fallback): Extracts only binaries from /opt.
#
# USAGE:
#   bash cache_ukrmol_docker_build.sh [IMAGE_NAME]
#
# EXAMPLE:
#   bash cache_ukrmol_docker_build.sh ukrmol-plus:latest
# ============================================================================ #

# 1. Configuration
IMAGE_NAME="${1:-ukrmol-plus:latest}"  # Default to 'ukrmol-plus:latest' if no arg provided
CACHE_DIR="./ukrmol-cache"
TEMP_CONTAINER_NAME="ukrmol-cache-extractor-tmp"

echo -e "\033[1;34m:: UKRmol+ Docker Cache Extractor ::\033[0m"
echo "Target Image: $IMAGE_NAME"
echo "Output Dir:   $CACHE_DIR"

# 2. Create Cache Directory
mkdir -p "$CACHE_DIR"

# 3. Create Temporary Container
# We create the container to access its filesystem, but we don't start it.
echo -e "\n\033[1;33mCreating temporary container from image...\033[0m"
# Remove any existing container with the same name just in case
docker rm "$TEMP_CONTAINER_NAME" > /dev/null 2>&1

if docker create --name "$TEMP_CONTAINER_NAME" "$IMAGE_NAME" > /dev/null 2>&1; then
    echo "Container created successfully."
else
    echo -e "\033[1;31mError: Could not create container. Does the image '$IMAGE_NAME' exist?\033[0m"
    exit 1
fi

# 4. Detect Build Stage (Builder vs Final)
# We check if the '/build' directory exists. 
# - If YES: It's a "Builder" stage image (contains downloads + source).
# - If NO:  It's a "Final" stage image (contains only compiled binaries).
echo -e "\n\033[1;33mDetecting image type...\033[0m"
# Attempt to list /build inside the container
IS_BUILDER=$(docker cp "$TEMP_CONTAINER_NAME":/build - 2>/dev/null | tar -t 2>/dev/null | head -1)

if [ -n "$IS_BUILDER" ]; then
    # --- STRATEGY A: BUILDER STAGE (Comprehensive Cache) ---
    echo ">> Detected 'Builder' stage layout."
    echo ">> Extracting ALL artifacts (Sources + Binaries)..."
    
    # Copy everything from /build to cache
    docker cp "$TEMP_CONTAINER_NAME":/build/. "$CACHE_DIR/"

else
    # --- STRATEGY B: FINAL STAGE (Binary Cache Only) ---
    echo ">> Detected 'Final' stage layout (Source code was discarded)."
    echo ">> Extracting compiled binaries only..."
    
    # The build script expects folders named: 'auxiliary-software' and 'ukrmolp-release'
    
    # 1. Compiler Toolchain
    echo "   - Extracting Compiler..."
    if docker cp "$TEMP_CONTAINER_NAME":/opt/compiler "$CACHE_DIR/auxiliary-software"; then
        echo "     [OK]"
    else
        echo "     [Warning] Could not find /opt/compiler"
    fi

    # 2. Compiled Libraries/Suite
    echo "   - Extracting UKRmol+ Suite..."
    if docker cp "$TEMP_CONTAINER_NAME":/opt/ukrmolp "$CACHE_DIR/ukrmolp-release"; then
        echo "     [OK]"
    else
        echo "     [Warning] Could not find /opt/ukrmolp"
    fi
fi

# 5. Cleanup
echo -e "\n\033[1;33mCleaning up...\033[0m"
# Remove git repos if they were copied (we want future builds to use your live git folders)
rm -rf "$CACHE_DIR/ukrmol-in-git"
rm -rf "$CACHE_DIR/ukrmol-out-git"
rm -rf "$CACHE_DIR/docker_build_script.sh"

# Remove the temp container
docker rm "$TEMP_CONTAINER_NAME" > /dev/null

echo -e "\n\033[1;32mSuccess! Cache populated in '$CACHE_DIR'.\033[0m"