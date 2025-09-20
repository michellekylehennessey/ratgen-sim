import click
import matplotlib.pyplot as plt
from .genetics import cross

@click.command()
@click.argument("sire")
@click.argument("dam")
@click.option("--litter", type=int, default=None, help="Simulate N pups.")
@click.option("--seed", type=int, default=None, help="Random seed for simulation.")
@click.option("--plot/--no-plot", default=True, help="Show phenotype probability bar chart.")
def main(sire, dam, litter, seed, plot):
    """
    Example:
      ratgen "A/a; P/p" "A/a; p/p" --litter 12
    """
    res = cross(sire, dam, litter_size=litter, seed=seed)
    click.echo("\nPhenotype probabilities:")
    for k, v in res["phenotypes"].items():
        click.echo(f"  {k:12s}  {v*100:5.1f}%")

    click.echo("\nGenotypes:")
    for k, v in res["table"].items():
        click.echo(f"  {k:16s}  {v*100:5.1f}%")

    if litter:
        click.echo(f"\nSimulated litter ({litter}): {', '.join(res['litter'])}")

    if plot:
        names = list(res["phenotypes"].keys())
        vals = [res["phenotypes"][n]*100 for n in names]
        plt.figure()
        plt.bar(names, vals)
        plt.ylabel("Probability (%)")
        plt.title(f"{sire} × {dam}")
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    main()
