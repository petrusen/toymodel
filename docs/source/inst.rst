Installation
------------

Here we detail a basic installation of the Toy Model (tested on Ubuntu 22):

.. code-block:: bash

   # (recommended) Create a dedicated environment
   python3 -m venv env4toymodel
   cd env4toymodel
   source bin/activate

   # Clone and install
   git clone https://github.com/petrusen/toymodel.git
   cd toymodel
   python3 -m pip install -e .  # dependencies already defined in setup.py
