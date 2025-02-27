# AABBtree

`AABBtree` is a library for the construction and manipulation of axis-aligned bounding box (AABB) trees. The library is designed to be used in the context of collision detection, and other applications that require the efficient computation of intersections between axis-aligned bounding boxes. It provides a simple and intuitive interface for the construction of AABB trees from a set of axis-aligned bounding boxes, as well as for the computation of intersections between AABB trees, individual AABBs, and points. The library is implemented in C++17 and is header-only, meaning that it can be easily integrated into existing C++ projects without the need for additional dependencies. This fast and efficient library is implemented in C++17 and is header-only, meaning that it can be easily integrated into existing C++ projects without the need for additional dependencies or compilation steps.

## Installation

### Quick and dirty

`AABBtree` is an header-only library that depends only on `CMake` (version >= 3.14).

### CMake

If you are using CMake, you can add the library as a subdirectory in your project.

```cmake
add_subdirectory(path/to/AABBtree)
target_link_libraries(your_target PRIVATE AABBtree::AABBtree)
```

You can use `FetchContent` to download the library from GitHub.

```cmake
include(FetchContent)

# Optionally specify a custom path to fetch content to
set(FETCHCONTENT_BASE_DIR "path/to/your/dependencies")
fetchcontent_declare(
  AABBtree
  GIT_REPOSITORY https://github.com/StoccoDavide/AABBtree.git
  GIT_TAG        main
)
fetchcontent_makeavailable(AABBtree)
target_link_libraries(your_target PRIVATE AABBtree::AABBtree)
```

If you already have `AABBtree` somewhere on your system, you can use `find_pacakge` directly.

```cmake
# Optionally specify a custom path to find content from
list(APPEND CMAKE_PREFIX_PATH "path/to/your/dependencies")
find_package(
  AABBtree
  ${YOUR_DESIRED_AABBTREE_VERSION}
  NO_MODULE
)

target_link_libraries(your_target PRIVATE AABBtree::AABBtree)
```

Since we are nice people, we also show you how to conditionally use `FetchContent` based if you already have the library or not.

```cmake
# Optionally specify a custom path to find content from
list(APPEND CMAKE_PREFIX_PATH "path/to/your/dependencies")
find_package(
  AABBtree
  ${YOUR_DESIRED_AABBTREE_VERSION}
  NO_MODULE
)

if(NOT TARGET AABBtree::AABBtree)
  include(FetchContent)

  # Optionally specify a custom path to fetch content to
  set(FETCHCONTENT_BASE_DIR "path/to/your/dependencies")
  fetchcontent_declare(
    AABBtree
    GIT_REPOSITORY https://github.com/StoccoDavide/AABBtree.git
    GIT_TAG        main
  )

  fetchcontent_makeavailable(AABBtree)
endif()

target_link_libraries(your_target PRIVATE AABBtree::AABBtree)
```

## Authors

- Davide Stocco <br>
  University of Trento <br>
  Department of Industrial Engineering <br>
  email: davide.stocco@unitn.it

- Enrico Bertolazzi <br>
  University of Trento <br>
  Department of Industrial Engineering <br>
  email: enrico.bertolazzi@unitn.it

Aka...

```
▗▄▄▄  ▄   ▄  ▐▌    ▗▞▀▜▌▄▄▄▄     ▐▌    ▗▄▄▖ ▗▞▀▚▖ ▄▄▄ ▄   ▄
▐▌  █ █   █  ▐▌    ▝▚▄▟▌█   █    ▐▌    ▐▌ ▐▌▐▛▀▀▘█    █   █
▐▌  █  ▀▄▀▗▞▀▜▌         █   █ ▗▞▀▜▌    ▐▛▀▚▖▝▚▄▄▖█     ▀▀▀█
▐▙▄▄▀     ▝▚▄▟▌               ▝▚▄▟▌    ▐▙▄▞▘          ▄   █
                                                       ▀▀▀
```

## License

The `AABBtree` project is distributed under the BSD 2-Clause License - see the [LICENSE](https://StoccoDavide.github.io/AABBtree/LICENSE) file for details.

Here's what the license entails:

1. Anyone can copy, modify and distribute this software.
2. You have to include the license and copyright notice with each and every distribution.
3. You can use this software privately.
4. You can use this software for commercial purposes.
5. This software is provided without warranty.
6. The software author or license can not be held liable for any damages inflicted by the software.
