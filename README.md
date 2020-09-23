# Variational Quantum Eigensolver Tutorial

## Installation

For best results:

```
conda create -n vqe
conda install jupyterlab qutip matplotlib
pip install qiskit[visualization]

```

and an account with IBM's quantum experience. Make sure to run: 

```
from qiskit import IBMQ
IBMQ.save_account('YOUR_TOKEN')
```

## Usage

```python vqe.py --samples 60 --n_shots 4096 --quantum/--classical/--backend name```

You can choose the number of samples of $\theta$, the number of shots per quantum circuit.
The default is to run locally. Use ```--classical``` to run on IBM's classical simulator, ```--quantum```
to run on ```ibmq_vigo```. 

For more details, check out the jupyter notebook!

