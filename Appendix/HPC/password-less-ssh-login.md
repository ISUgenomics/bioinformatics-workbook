# Setup password-less login for HPC

To login automatically from your machine to the remote host, you can save the private/public key pair in both machines. This way, you don't have to enter password each time you login.

## Step 1: Create public and private keys (local computer) </h3>

```
ssh-keygen
```

## Step 2: Copy the public key to remote-host </h3>

```
ssh-copy-id -i ~/.ssh/id_rsa.pub userid@lightning3.its.iastate.edu
```

## Step 3: Login  </h3>

```
ssh userid@lightning3.its.iastate.edu
# Ensure file permissions for ~/.ssh/.id_rsa (local) and ~/.ssh/authorized_keys (remote) are such that it is only readable by you!
```
