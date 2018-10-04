# Setup password-less login for HPC

Most HPC resources should now require double authentication which will make this type of password-less login not possible.  For those that don't have double authentication, this will work.

To login automatically from your machine to the remote host, you can save the private/public key pair in both machines. This way, you don't have to enter password each time you login.

## Step 1: Create public and private keys (local computer) </h3>

```
ssh-keygen
```

## Step 2: Copy the public key to remote-host </h3>

```
ssh-copy-id -i ~/.ssh/id_rsa.pub userid@MACHINENAME
```

## Step 3: Login  </h3>

```
ssh userid@MACHINENAME
# Ensure file permissions for ~/.ssh/.id_rsa (local) and ~/.ssh/authorized_keys (remote) are such that it is only readable by you!
```
