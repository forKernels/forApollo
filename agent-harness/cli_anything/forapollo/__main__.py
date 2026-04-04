"""Entry point for cli-anything-forapollo."""

from .forapollo_cli import cli

def main():
    cli(standalone_mode=False)

if __name__ == "__main__":
    cli()
