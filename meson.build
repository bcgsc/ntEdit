project(
    'ntEdit',
    'cpp',
    version: '2.0.0',
    license: 'GPL-3',
    default_options: [
        'cpp_std=c++17',
        'werror=true',
        'warning_level=3',
        'optimization=3',
        'b_lto=true',
    ]
)

compiler = meson.get_compiler('cpp')

dependencies = [
    dependency('openmp'),
    dependency('zlib'),
    compiler.find_library('btllib'),
    compiler.find_library('nthash'),
]

executable(
    'ntedit',
    'ntedit.cpp',
    dependencies: dependencies,
    include_directories: 'lib',
    install: true,
)
