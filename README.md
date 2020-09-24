# Variational Quantum Eigensolver Tutorial

This is my submission (Task 4) for the QOSF Mentorship program. 

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

You can choose the number of samples of &Theta; as well as the number of shots per quantum circuit.
The default is to run locally. Use ```--classical``` to run on IBM's classical simulator, ```--quantum```
to run on ```ibmq_vigo```. Or choose your own backend.

For more details, check out the jupyter notebook!

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/heyredhat/vqe_tutorial/master?filepath=VQE_Tutorial.ipynb)

