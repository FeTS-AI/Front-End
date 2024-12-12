import argparse
from FigureGenerator.screenshot_maker import figure_generator

def setup_argparser():
    parser = argparse.ArgumentParser("Figure Generator Wrapper")
    parser.add_argument("-i", type=str)
    parser.add_argument("-ylabels", type=str)
    parser.add_argument("-output", type=str)
    parser.add_argument("-masks", type=str, nargs='?', const=None, default=None)
    return parser

def main():
    args = setup_argparser().parse_args()
    images = args.i
    ylabels = args.ylabels
    output_filename = args.output
    input_mask = args.masks

    figure_generator(
        input_images=images,
        ylabels=ylabels,
        output=output_filename,
        input_mask=input_mask,
        flip_sagittal=True,
        flip_coronal=True,
    )

if __name__ == "__main__":
    main()