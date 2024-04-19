project(
    'ntEdit',
    'cpp',
    version: '2.0.1',
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
    compiler.find_library('boost_iostreams'),
]

include_dirs = [
    include_directories('src/vendor/argparse/include'),
]

link_args = [ '-lbtllib' ]

executable(
    'ntedit',
    'ntedit.cpp',
    dependencies: dependencies,
    include_directories: 'lib',
    link_args: link_args,
    install: true,
)

src_files = [ 'ntedit_make_genome_bf.cpp' ]
# Source files and scripts
subdir('src')

scripts = ['run-ntedit', 'ntedit_run_pipeline.smk', 'ntedit_cross_reference_vcf.py']
install_data(scripts, install_dir : 'bin', install_mode : 'rwxr-xr--')
