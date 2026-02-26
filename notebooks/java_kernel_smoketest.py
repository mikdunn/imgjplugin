from __future__ import annotations

import time

from jupyter_client import KernelManager


def main() -> int:
    km = KernelManager(kernel_name="java")
    km.start_kernel()
    kc = km.client()
    kc.start_channels()

    try:
        kc.execute('System.out.println("hello from java kernel"); int a=1; int b=2; System.out.println(a+b);')

        got_idle = False
        streams: list[str] = []
        results: list[str] = []
        errors: list[str] = []

        t0 = time.time()
        while time.time() - t0 < 30:
            try:
                msg = kc.get_iopub_msg(timeout=1)
            except Exception:
                continue

            mtype = msg.get("header", {}).get("msg_type")
            content = msg.get("content", {})

            if mtype == "stream":
                streams.append(str(content.get("text", "")))
            elif mtype == "execute_result":
                data = content.get("data", {})
                results.append(str(data.get("text/plain", "")))
            elif mtype == "error":
                errors.append(str(content.get("evalue", "")))
            elif mtype == "status" and content.get("execution_state") == "idle":
                got_idle = True
                break

        print(f"idle: {got_idle}")
        for s in streams:
            s = s.strip()
            if s:
                print(f"stream: {s}")
        for r in results:
            r = r.strip()
            if r:
                print(f"result: {r}")
        for e in errors:
            e = e.strip()
            if e:
                print(f"error: {e}")

        return 0 if got_idle and not errors else 1
    finally:
        try:
            kc.stop_channels()
        except Exception:
            pass
        try:
            km.shutdown_kernel(now=True)
        except Exception:
            pass


if __name__ == "__main__":
    raise SystemExit(main())
