/* Copyright 2020
 *
 * Zdenek Masin with contributions from others (see the UK-AMOR website)
 *
 * This file is part of GBTOlib.
 *
 *     GBTOlib is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     GBTOlib is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
 *     <https://www.gnu.org/licenses/>.
 */

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32

#include <objbase.h>
#include <windows.h>

#else

#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

#include <sys/mman.h>
#include <sys/stat.h>

#endif

/// \brief   Map data segment from file
/// \authors J Benda
/// \date    2020 - 2025
///
/// Assigns a virtual memory range to a contiguous data segment within a disk file, which can be used as any other memory
/// buffer. The file need not exist; if `create_new` is set to 1, it will be created first, sufficiently large to hold the
/// requested memory range. The file name can be left partially unspecified (by setting `unique_name` to 1) to allow the
/// operating system to generate a unique file path with the prefix given by `filename`. Such automatically created files
/// will be deleted once the mapping is canceled (or the program ends).
///
/// \param[in]  filename    Path to the file to map or a prefix for automatically generated name.
/// \param[in]  offset      Byte position in file to begin the mapping.
/// \param[in]  length      Byte length of the mapped window (buffer).
/// \param[out] addr        Internal address that needs to be passed to `unmap_file_window` to cancel the mapping.
/// \param[out] ptr         Pointer to memory where the mapped data begin.
/// \param[in]  read_write  If set to 1, allow write access to the mapped memory.
/// \param[in]  unique_name If set to 1, interpret `filename` only as a prefix, and supplement it with a unique ending.
/// \param[in]  create_new  If set to 1, create a new file instead of mapping an existing one.
///
void map_file_window
(
    const char* filename,
    size_t offset,
    size_t length,
    void** addr,
    void** ptr,
    int read_write,
    int unique_name,
    int create_new
)
{
    const char* debug = getenv("GBL_DEBUG_MMAP");

    int name_length = strlen(filename);
    int extra_length = 256;

#ifdef _WIN32
    char* fullname = (char*)malloc(MAX_PATH);
#else
    char* fullname = (char*)malloc(name_length + extra_length + 1);
#endif

    strcpy(fullname, filename);

    if (debug) printf("[map_file_window]: filename = \"%s\"\n", name_length > 0 ? filename : "(anonymous)");
    if (debug) printf("[map_file_window]: offset = %ld\n", offset);
    if (debug) printf("[map_file_window]: length = %ld\n", length);

    // Get virtual memory page size
    //  - Data segments mapped from file to memory must begin on page boundary, i.e. integer multiples
    //    of the page size. This is typically 4 kiB.

#ifdef _WIN32
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    long pagesize = si.dwAllocationGranularity;
#else
    long pagesize = sysconf(_SC_PAGE_SIZE);
#endif
    if (debug) printf("[map_file_window]: page size = %ld\n", pagesize);

    // Create unique file name or template
    //  - If unique name is requested, we use different strategies for creation of scratch files in the provided path.
    //    In Windows we generate a unique file name based on the prefix path.
    //    In Linux we prepare template for the library call `mkstemp`.

    if (unique_name)
    {
#ifdef _WIN32
        // split filename to directory and file name prefix (search separator from end)
        const char* prefix = filename + strlen(filename);
        while (prefix > filename && *(prefix - 1) != '/' && *(prefix - 1) != '\\')
            prefix--;
        char* directory = strdup(filename);
        directory[prefix - filename] = 0;
        for (char* back = directory + (prefix - filename) - 1; back >= directory && (*back == '/' || *back == '\\'); back--)
            *back = 0;
        // get unique file name based on the input
        GetTempFileName(directory, prefix, 0, fullname);
        // clean up
        free(directory);
#else
        strcat(fullname, "XXXXXX");
#endif
    }

    // Open the backing file
    //  - Open the disk file (possibly creating a new one) that will be mapped from.

#ifdef _WIN32
    LARGE_INTEGER li;
    HANDLE hFile = INVALID_HANDLE_VALUE;
    if (name_length > 0)
    {
        hFile = CreateFile
        (
            fullname,
            (read_write ? GENERIC_WRITE : 0 ) | GENERIC_READ,
            FILE_SHARE_READ | FILE_SHARE_WRITE,
            NULL,
            create_new ? CREATE_ALWAYS : OPEN_EXISTING,
            create_new ? (unique_name ? FILE_FLAG_DELETE_ON_CLOSE : 0) | FILE_ATTRIBUTE_TEMPORARY : FILE_ATTRIBUTE_NORMAL,
            NULL
        );
        if (hFile == INVALID_HANDLE_VALUE)
        {
            printf("Failed to %s file \"%s\", error %d.\n", create_new ? "create" : "open", fullname, GetLastError());
            exit(1);
        }
    }
#else
    int fd = -1;
    if (name_length > 0 || (create_new && unique_name))
    {
        if (create_new)
        {
            fd = unique_name ? mkstemp(fullname) : open(filename, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
        }
        else
        {
            fd = open(filename, read_write ? O_RDWR : O_RDONLY);
        }

        if (fd < 0)
        {
            printf("Failed to open file \"%s\", error %d.\n", fullname, errno);
            exit(1);
        }
    }
#endif

    if (debug) printf("[map_file_window]: fullname = \"%s\"\n", fullname);

    // Resize the underlying new file
    //  - Here we reserve space for the mapping by resizing any newly created scratch file. In both systems this results
    //    in a sparse file which has a nominal size but does not really occupy any disk space. But this is all that is needed
    //    for the mapping routines.

    if (create_new)
    {
#ifdef _WIN32
        li.QuadPart = offset + length;
        if (!SetFilePointerEx(hFile, li, NULL, FILE_BEGIN))
        {
            printf("Failed to set pointer to offset %ld in file \"%s\", error %d.\n", offset + length, fullname, GetLastError());
            exit(1);
        }
        if (!SetEndOfFile(hFile))
        {
            printf("Failed to resize file \"%s\" to offset %ld, length %ld, error %d.\n", fullname, offset, length, GetLastError());
            exit(1);
        }
#else
        if (fd >= 0 && ftruncate(fd, offset + length) != 0)
        {
            printf("Failed to resize file \"%s\" to offset %ld, length %ld, error %d.\n", fullname, offset, length, errno);
            exit(1);
        }
#endif
    }

    // Align offset to page boundary
    //  - Find position in file less than or equal to the provided offset, which is aligned to the page boundary.
    //    This will be the actual beginning of our mapping.

    size_t pa_offset = offset & ~(pagesize - 1);
    if (debug) printf("[map_file_window]: pa_offset = %ld\n", pa_offset);

    // Create the mapping
    //  - Assign virtual memory range to the segment from `pa_offset` to `offset` + `length`.

#ifdef _WIN32
    li.QuadPart = offset + length - pa_offset;
    HANDLE hFileMapping = name_length == 0 ?
        CreateFileMapping(hFile, NULL, read_write ? PAGE_READWRITE : PAGE_READONLY, li.HighPart, li.LowPart, NULL) :
        CreateFileMapping(hFile, NULL, read_write ? PAGE_READWRITE : PAGE_READONLY, 0, 0, NULL);
    if (hFileMapping == NULL)
    {
        printf("Failed to create mapping for file \"%s\", error %d.\n", fullname, GetLastError());
        exit(1);
    }
    li.QuadPart = pa_offset;
    *addr = MapViewOfFile(hFileMapping, read_write ? FILE_MAP_WRITE : FILE_MAP_READ, li.HighPart, li.LowPart, length + offset - pa_offset);
    if (*addr == NULL)
    {
        printf("Failed to map memory from file \"%s\", error %d.\n", filename, GetLastError());
        exit(1);
    }
#else
    *addr = mmap
    (
        NULL,
        length + offset - pa_offset,
        read_write ? PROT_READ | PROT_WRITE : PROT_READ,
        fd >= 0 ? MAP_SHARED : MAP_PRIVATE | MAP_ANONYMOUS,
        fd,
        pa_offset
    );
    if (*addr == MAP_FAILED)
    {
        printf("Failed to map file \"%s\" of size %ld from offset %ld: %s\n", fullname, length, offset, strerror(errno));
        exit(1);
    }

#endif

    // Point to the requested offset
    //  - Set up a pointer that actually points to the position `offset` in the file. This is generally not identical
    //    to the beginning of the mapping. But for the user it will be the beginning. In case of a newly create file with
    //    zero offset, the mapping always starts at the end of file, which is already page-aligned by default.

    *ptr = *addr + offset - pa_offset;

    // Close the backing file
    //  - This closes the underlying file, but keeps its contents mapped to memory. In case of a scratch file (`unique_name` = 1),
    //    we even delete the file altogether (by `unlink` in Linux, or via `FILE_FLAG_DELETE_ON_CLOSE` in Windows). This will
    //    make sure that we do not need to clean it up afterwards when the mapping is finalized (or the program crashes).

#ifdef _WIN32
    CloseHandle(hFileMapping);
    CloseHandle(hFile);
#else
    if (unique_name)
    {
        unlink(fullname);
    }
    close(fd);
#endif

    if (debug && !create_new)
    {
        printf("[map_file_window]: first double = %g\n", *(double*)*ptr);
        printf("[map_file_window]: first long = %ld\n", *(long*)*ptr);
    }
}


/// \brief   Cancel file mapping
/// \authors J Benda
/// \date    2020 - 2021
///
/// Release memory mapping, allowing e.g. final deletion of the file.
///
/// \param[in] addr    Internal address obtained from `map_file_window`.
/// \param[in] length  Length of the mapped data segment.
///
void unmap_file_window (void* addr, size_t length)
{
#ifdef _WIN32
    UnmapViewOfFile(addr);
#else
    munmap(addr, length);
#endif
}
