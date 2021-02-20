.. _Continuous_Integration_process:

Building GEOSX using Docker with precompiled TPL binaries
=========================================================

For development purposes, you may want to use the publicly available docker images or the OSX tarball instead of compiling them yourself.
While this is surely possible, please note that *this is not supported by the GEOSX team that reserves the right to modify its workflow or delete elements on which you may have build your own workflow*.

There are multiple options to use the exposed docker images.

- A lot of IDE now provide remote development modes (e.g. `CLion <https://www.jetbrains.com/help/clion/remote-projects-support.html>`_, `VS Code <https://code.visualstudio.com/docs/remote/remote-overview>`_, `Eclipse Che <https://www.eclipse.org/che/>`_ and surely others).
  Depending on your choice, please read their documentation carefully so you can add their own requirements on top the TPL images that are already available.
  Feel free to share them :)
- Another option is to develop directly inside the container (*i.e.* not remotely).
  Install your favorite development inside the image (be mindful of X display issues), connect to the running container and start hacking!

Make sure to build GEOSX with the environment variables provided in the :ref:`Docker_images_contract`.

Please be aware of how to retrieve back your code: you may want to bind mount volumes (``-v``/``--volume=`` options of `docker run <https://docs.docker.com/engine/reference/run/>`_).
Think about opening ports too.

Getting a Docker container
--------------------------

To obtain a Docker image follow the three steps below:

1. `Install <https://docs.docker.com/get-docker/>`_ Docker on your machine
2. `Select <https://hub.docker.com/u/geosx/>`_ your favorite image
3. Use the ``docker run Image[:tag]`` command specifying the image to derive the container from

For example,

.. code-block:: sh

    docker login
    docker run -ti docker run -ti geosx/ubuntu18.04-gcc8:140-613

will derive the container for a ubuntu18.04 OS using the gcc-8 compiler suite, where the tag `:140-613` specifies the Travis pull request and build number, respectively.
For additional details see the `Docker run reference page <https://docs.docker.com/engine/reference/run/>`_.

Follow the instructions provided in the Quick Start Guide page to download, configure and compile GEOSX. 
