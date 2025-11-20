from pystencilssfg import SfgComposer
from pystencils import Target
from sweepgen import get_build_config
import argparse


def configure(sfg: SfgComposer, add_defines: bool = False):
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--platform", dest="platform", default="cpu", type=str)

    args = parser.parse_args(sfg.context.argv)

    build_config = get_build_config(sfg)

    match args.platform:
        case "cpu":
            build_config.override.target = Target.CurrentCPU

            if add_defines:
                sfg.code("#define CouetteFlow_CPU_BUILD true")
        case "hip":
            build_config.override.target = Target.HIP

            if add_defines:
                sfg.code("#define CouetteFlow_GPU_BUILD true")
        case "cuda":
            build_config.override.target = Target.CUDA

            if add_defines:
                sfg.code("#define CouetteFlow_GPU_BUILD true")
        case _:
            raise ValueError(f"Unexpected target id: {args.target}")
