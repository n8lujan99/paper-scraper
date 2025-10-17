import os, smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from dotenv import load_dotenv

# load local .env if present
load_dotenv()

def send_email(cfg, subject, text_body, html_body):
    """
    send an email using credentials from environment variables
    or config.yaml (via cfg["output"]["email"]).

    the YAML config should define:
      output:
        email:
          from_addr: ...
          to_addrs: [...]
          username: ...
          password_env: "EMAIL_PASS"   # environment variable name
          smtp_host: "smtp.gmail.com"
          smtp_port: 587
          use_starttls: true
    """
    em = cfg["output"]["email"]

    # build MIME message
    msg = MIMEMultipart("alternative")
    msg["subject"] = subject
    msg["from"] = em["from_addr"]
    msg["to"] = ", ".join(em["to_addrs"])
    msg.attach(MIMEText(text_body, "plain", "utf-8"))
    msg.attach(MIMEText(html_body, "html", "utf-8"))

    # load password securely
    password = os.getenv(em.get("password_env", "EMAIL_PASS"), "")
    if not password:
        raise RuntimeError(
            f"missing password in environment variable: {em.get('password_env', 'EMAIL_PASS')}"
        )

    # connect and send
    with smtplib.SMTP(em["smtp_host"], em["smtp_port"]) as s:
        if em.get("use_starttls", True):
            s.starttls()
        s.login(em["username"], password)
        s.sendmail(em["from_addr"], em["to_addrs"], msg.as_string())
        print(f"email sent to {', '.join(em['to_addrs'])}")
