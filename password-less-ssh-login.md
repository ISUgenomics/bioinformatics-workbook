To login automatically from your machine to the remote host, you can save the private/public key pair in both machines. This way, you don't have to enter password each time you login.

<h3> Step 1: Create public and private keys (local computer) </h3>

<pre>
ssh-keygen
</pre>

<h3> Step 2: Copy the public key to remote-host </h3>

<pre>
ssh-copy-id -i ~/.ssh/id_rsa.pub userid@lightning3.its.iastate.edu
</pre>

<h3> Step 3: Login  </h3>

<pre>
ssh userid@lightning3.its.iastate.edu
# Ensure file permissions for ~/.ssh/.id_rsa (local) and ~/.ssh/authorized_keys (remote) are such that it is only readable by you! 
</pre>