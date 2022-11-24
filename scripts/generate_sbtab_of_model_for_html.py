#!/usr/bin/env python3
"""
Import model and generate SBtab of model structure.
Arguments: 
    model_dir: (relative) path to RBA-model
    --output-dir: (relative) path to storage location of results.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import argparse
import sys

# package imports
from rbatools.rba_session import SessionRBA

def main():
    parser = argparse.ArgumentParser(description='Generate SBtabs (with links) for RBA model and also for sample outputs')
    parser.add_argument('model_dir', metavar='model-dir', type=str,
                        help='Directory of RBA-model')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Directory to save the SBtab file')

    args = parser.parse_args()
    Session=SessionRBA(xml_dir=args.model_dir,lp_solver="swiglpk")

    if args.output_dir is not None:
        output_dir=args.output_dir+"/"
    else:
        output_dir=""
        
    Session.ModelStructure.export_sbtab(filename=output_dir+"RBAmodel_SBtab_with_links",add_links=True)
    print("Model SBtab written to: {}".format(output_dir+"RBAmodel_SBtab_with_links.tsv"))

if __name__ == '__main__':
    main()
